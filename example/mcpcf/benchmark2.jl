
# benchmark with mc results
# use BSeq to compute and compare

using ElectronGas
using ElectronGas.Parameter
using ElectronGas.Convention
using ElectronGas.LegendreInteraction
using ElectronGas.Interaction
using ElectronGas.Parameters
using ElectronGas.GreenFunc
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids
using ElectronGas.SelfEnergy
using ElectronGas.JLD2
using ElectronGas.BSeq: initFR, calcF!, calcR!, G02wrapped
using Test

function oneloop_solver(param,
    G2::GreenFunc.MeshArray, kernel, kernel_ins,
    qgrids::Vector{CompositeGrid.Composite}, Euv;
    source::Union{Nothing,GreenFunc.MeshArray}=nothing,
    source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1], G2.mesh[2]; dtype=Float64, data=ones(1, G2.mesh[2].size)),
    verbose=false)

    if verbose
        println("atol=$atol,rtol=$rtol")
    end

    @unpack dim, kF = param
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)
    if !(source isa Nothing)
        @assert source.mesh[1] isa MeshGrids.ImTime "ImTime is expect for the dim = 1 source."
        for ni in eachindex(F_freq.mesh[1])
            source[ni, :] .*= kgrid.grid
        end
    end
    source_ins[1, :] .*= kgrid.grid
    source0 = source_ins[1, kF_label]

    # Initalize F and R
    F_freq, R_imt, R_ins = initFR(Euv, G2.mesh[1].representation.rtol, kgrid, param)
    R_ins.data .= 1.0
    R_imt.data .= 0.0

    R0, R0_sum = 1.0, 0.0
    dR0, dR0_sum = zeros(Float64, (kgrid.size)), zeros(Float64, (kgrid.size))
    R_sum = zeros(Float64, (R_imt.mesh[1].representation.size, kgrid.size))

    calcF!(F_freq, R_imt, R_ins, G2)

    if dim == 3
        calcR!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
    elseif dim == 2
        calcR_2d!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
    else
        error("Not implemented for $dim dimensions.")
    end

    R_kF = real(dlr_to_imfreq(to_dlr(R_imt), [0])[1, kF_label] + R_ins[1, kF_label])
    # split the low-energy part R0=R(ω₀,kF) and the remaining instant part dR0 for iterative calcualtions 
    R0_sum = source0 + R_kF
    dR0_sum = view(R_ins, 1, :) + view(source_ins, 1, :) .- (source0 + R_kF)
    R0 = R0_sum
    dR0 = dR0_sum
    R_ins[1, :] = dR0 .+ R0

    # iterative calculation of the dynamical part 
    R_sum = view(R_imt, :, :)
    R_imt[:, :] = R_sum
    @debug "R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))"
    # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

    for ni in eachindex(F_freq.mesh[1])
        F_freq[ni, :] ./= kgrid.grid
        R_imt[ni, :] ./= kgrid.grid
    end
    R_ins[1, :] ./= kgrid.grid

    return F_freq, R_imt, R_ins
end

function oneloop(param, channel::Int;
    Euv=100 * param.EF, rtol=1e-10, atol=1e-10,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8,
    Vph::Union{Function,Nothing}=nothing, sigmatype=:none, int_type=:rpa,
    α=0.8, verbose=false, Ntherm=30, Nmax=10000,
    issave=false, uid=1, dir="./")
    @unpack dim, rs, β, kF = param
    if verbose
        println("atol=$atol,rtol=$rtol")
    end
    # prepare Legendre decomposed effective interaction
    if dim == 3
        if channel == 0
            # nmax = G2.mesh[1].grid[end]
            # delta_correction = Δω_correction(G2, wsph, a2f_iso, nmax)
            @time W = LegendreInteraction.DCKernel0(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph)
        else
            @time W = LegendreInteraction.DCKernel_old(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph)
        end
    elseif dim == 2
        @time W = LegendreInteraction.DCKernel_2d(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
            minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph)
    else
        error("No support for $dim dimension!")
    end

    fdlr = Lehmann.DLRGrid(Euv, β, rtol, true, :pha)
    bdlr = W.dlrGrid
    kgrid = W.kgrid
    qgrids = W.qgrids

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_ins = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3))

    kF_label = locate(kgrid, kF)
    qF_label = locate(qgrids[kF_label], kF)
    println("static kernel at (kF, kF):$(kernel_ins[kF_label, qF_label])")

    # prepare G2 as sigmatype
    if sigmatype == :none
        G2 = G02wrapped(Euv, rtol, kgrid, param)
    elseif sigmatype == :g0w0
        @unpack me, β, μ = param
        wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
        Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        # self energy should be converted to proper frequency grid
        Σ_dlr = Σ |> to_dlr
        Σ_wn = dlr_to_imfreq(Σ_dlr, wn_mesh)
        G2 = G2wrapped(Σ_wn, Σ_ins, param)
    end

    # calculate F, R by Bethe-Slapter iteration.
    F_freq, R_imt, R_ins = oneloop_solver(param, G2, kernel, kernel_ins, qgrids, Euv; verbose=verbose)

    R_freq = R_imt |> to_dlr |> to_imfreq
    if issave
        fname = "oneloopdata_$(uid).jld2"
        jldopen(dir * fname, "w") do file
            file["param"] = param
            file["F_freq"] = F_freq
            file["R_ins"] = R_ins
            file["R_freq"] = R_freq
        end
    end
    # println(view(R_freq, :, kF_label))
    return R_freq, F_freq, R_ins
end

function R0(ri, rw, param)
    kF = param.kF
    kgrid = rw.mesh[2]
    ikF = searchsortedfirst(kgrid, kF)
    # return 1.0 .+ ri[ikF] .+ rw[:, ikF]
    return ri[ikF] .+ rw[:, ikF]
end

param = Parameter.defaultUnit(0.1, 3.0)
Rw, Fw, Ri = oneloop(param, 0;
    int_type=:rpa)

println(Ri.mesh[2])
println(Ri.data)
println(R0(Ri, Rw, param))