"""
Bethe-Slapter-type equation solver and the application to Cooper-pair linear response approach. 
"""
module BSeq_resum

using ..Parameter, ..Convention, ..LegendreInteraction, ..Interaction
using ..Parameters, ..GreenFunc, ..Lehmann, ..CompositeGrids
using ..SelfEnergy
using ..JLD2
using ..BSeq

const freq_sep = 0.01

function initFR_resum(Euv, rtol, sgrid, param)
    @unpack β, me, μ = param
    Ω_c = freq_sep

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    R_freq = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
    F_freq = similar(R_freq)
    for ind in eachindex(R_freq)
        ni, ki = ind[1], ind[2]
        ωn, k = wn_mesh[ni], sgrid[ki]
        e = k^2 / 2 / me - μ
        R_freq[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
    end
    R_imt = real(R_freq |> to_dlr |> to_imtime)
    R_ins = GreenFunc.MeshArray([1], sgrid; dtype=Float64, data=zeros(1, sgrid.size))
    F_fs = GreenFunc.MeshArray(wn_mesh, [1,]; dtype=Float64, data=zeros(length(wn_mesh), 1))
    return F_freq, F_fs, R_imt, R_ins
end

function initB_resum(W, Euv, rtol, param)
    @unpack dim, kF = param

    sgrid = W.kgrid
    kgrid = sgrid
    qgrids = W.qgrids
    bdlr = W.dlrGrid
    iqFs = [locate(qgrids[ki], kF) for ki in 1:kgrid.size]


    kernel_ins = W.kernel_bare
    kernel_freq = W.kernel

    @unpack β, me, μ = param
    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    B = GreenFunc.MeshArray(wn_mesh, wn_mesh, sgrid; dtype=Float64)
    W0 = similar(B)

    for ind in eachindex(W0)
        in1, in2, ik = ind[1], ind[2], ind[3]
        w1, w2, k = wn_mesh[in1], wn_mesh[in2], sgrid[ik]
        n1, n2 = wn_mesh.grid[in1], wn_mesh.grid[in2]
        # V = kernel_ins[ik, iqFs[ik]]
        Ww = view(kernel_freq, ik, iqFs[ik], :)
        Wdlr = Lehmann.matfreq2dlr(bdlr, Ww)
        W12 = real(Lehmann.dlr2matfreq(bdlr, Wdlr, [n1 - n2,]))
        # println((n1, n2))
        # println((V, W12[1]))
        # W0[ind] = V + W12[1]
        W0[ind] = W12[1]
    end
    B.data .= W0.data
    W0dlr = GreenFunc.to_dlr(W0; dim=2)
    W0 = GreenFunc.dlr_to_imtime(W0dlr; dim=2)
    return B, W0
end

function calcF_resum!(F::GreenFunc.MeshArray, F_fs::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray, Πs::GreenFunc.MeshArray; ikF)
    R_freq = R |> to_dlr |> to_imfreq
    # algebraic in frequency space
    for ind in eachindex(F)
        F[ind] = real(R_freq[ind] + R_ins[1, ind[2]]) * G2[ind]
    end

    for ind in eachindex(F_fs)
        F_fs[ind] = -real(R_freq[ind[1], ikF]
                          +
                          R_ins[1, ikF]) * Πs[ind]
    end
    # for ind in 1:length(F_fs.mesh[1])
    #     F_fs[ind, 1] = -real(R_freq[ind, ikF]
    #                          +
    #                          R_ins[1, ikF]) * Πs[ind, 1]
    # end
end

function calcF_freq_resum!(F::GreenFunc.MeshArray, F_fs::GreenFunc.MeshArray, R::GreenFunc.MeshArray, G2::GreenFunc.MeshArray, Πs::GreenFunc.MeshArray; ikF)
    R_freq = R # assume R_ins is included
    # algebraic in frequency space
    for ind in eachindex(F)
        F[ind] = real(R_freq[ind]) * G2[ind]
    end

    for ind in eachindex(F_fs)
        F_fs[ind] = -real(R_freq[ind[1], ikF]) * Πs[ind]
    end
end

function calcR_resum!(F::GreenFunc.MeshArray, F_fs::GreenFunc.MeshArray,
    R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
    source::Union{Nothing,GreenFunc.MeshArray},
    kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite}; iqFs)
    kgrid = F.mesh[2]
    F_dlr = F |> to_dlr
    # switch to τ space
    F_imt = real(F_dlr |> to_imtime)
    F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)

    Ffs_dlr = F_fs |> to_dlr
    # switch to τ space
    Ffs_imt = real(Ffs_dlr |> to_imtime)
    Ffs_ins = real(dlr_to_imtime(Ffs_dlr, [F_fs.mesh[1].representation.β,])) * (-1)

    for ind in eachindex(R)
        # for each τ, k, integrate over q
        τi, ki = ind[1], ind[2]
        # interpolate F to q grid of given k
        Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
        integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq
        R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        R[τi, ki] += kernel[ki, iqFs[ki], τi] * Ffs_imt[τi, 1] ./ (-4 * π * π)
        if τi == 1
            # same for instant part
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
            integrand = view(kernel_ins, ki, 1:qgrids[ki].size) .* Fq
            R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            R_ins[1, ki] += kernel_ins[ki, iqFs[ki]] * Ffs_ins[1, 1] ./ (-4 * π * π)
        end
    end
    !(source isa Nothing) && (R += source)
end

function Rt2Rw!(Rw::GreenFunc.MeshArray, Rt::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray)
    R_freq = Rt |> to_dlr |> to_imfreq
    for ind in eachindex(Rw)
        Rw[ind] = real(R_freq[ind] + R_ins[1, ind[2]])
    end
end

function Πs0wrapped(Euv, rtol, param; ω_c=0.02param.EF)
    @unpack me, β, μ, kF, EF = param

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    green = GreenFunc.MeshArray(wn_mesh, [1,]; dtype=Float64)
    for ind in eachindex(green)
        ni, ki = ind[1], ind[2]
        ωn = wn_mesh[ni]
        # ω1 = π * param.T
        # green[ind] = π * me / (kF) / abs(ωn) / (1 + (abs(ωn) / (ω_c))^2) * (1 + (ω1 / ω_c)^2)
        green[ind] = 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    end
    return green
end

function BSeq_solver_resum(param, G2::GreenFunc.MeshArray, Πs::GreenFunc.MeshArray,
    kernel, kernel_ins,
    qgrids::Vector{CompositeGrid.Composite},
    Euv; Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8,
    source::Union{Nothing,GreenFunc.MeshArray}=nothing,
    source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1], G2.mesh[2];
        dtype=Float64, data=ones(1, G2.mesh[2].size)),
    verbose=false, Ncheck=5, Nmax=10000)

    if verbose
        println("atol=$atol,rtol=$rtol")
    end

    @unpack dim, kF = param
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)
    ikF = kF_label
    iqFs = [locate(qgrids[ki], kF) for ki in 1:kgrid.size]
    println("kF=$kF")
    # println("ikF=$ikF")
    # println("iqFs=$iqFs")

    if !(source isa Nothing)
        @assert source.mesh[1] isa MeshGrids.ImTime "ImTime is expect for the dim = 1 source."
        for ni in eachindex(source.mesh[1])
            source[ni, :] .*= kgrid.grid
        end
    end
    source_ins[1, :] .*= kgrid.grid
    source0 = source_ins[1, kF_label]

    # Initalize F and R
    F_freq, F_fs, R_imt, R_ins = initFR_resum(Euv, G2.mesh[1].representation.rtol, kgrid, param)

    # R0, R0_sum = 1.0, 0.0
    R0, R0_sum = source0, 0.0
    dR0, dR0_sum = zeros(Float64, (kgrid.size)), zeros(Float64, (kgrid.size))
    R_sum = zeros(Float64, (R_imt.mesh[1].representation.size, kgrid.size))

    lamu, lamu0 = 0.0, 1.0
    n = 0
    # self-consistent iteration with mixing α
    # Note!: all quantites about R, F, and source in the `while` loop are multiplied by momentum k.
    while (true)
        n = n + 1
        # switch between imtime and imfreq to avoid convolution
        # dlr Fourier transform is much faster than brutal force convolution

        # calculation from imtime R to imfreq F
        calcF_resum!(F_freq, F_fs, R_imt, R_ins, G2, Πs; ikF=ikF)

        # calculation from imfreq F to imtime R
        if dim == 3
            calcR_resum!(F_freq, F_fs, R_imt, R_ins, source, kernel, kernel_ins, qgrids; iqFs=iqFs)
            # elseif dim == 2
            #     calcR_2d!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
        else
            error("Not implemented for $dim dimensions.")
        end

        R_kF = real(dlr_to_imfreq(to_dlr(R_imt), [0])[1, kF_label] + R_ins[1, kF_label])
        # split the low-energy part R0=R(ω₀,kF) and the remaining instant part dR0 for iterative calcualtions 
        R0_sum = source0 + R_kF + R0_sum * α
        dR0_sum = view(R_ins, 1, :) + view(source_ins, 1, :) .- (source0 + R_kF) + dR0_sum .* α
        R0 = R0_sum * (1 - α)
        dR0 = dR0_sum .* (1 - α)
        R_ins[1, :] = dR0 .+ R0

        # iterative calculation of the dynamical part 
        R_sum = view(R_imt, :, :) + R_sum .* α
        R_imt[:, :] = R_sum .* (1 - α)
        @debug "R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))"
        # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

        # record lamu=1/R0 if iterative step n > Ntherm
        if n > Ntherm && (n % Ncheck == 1)
            lamu = -1 / (1 + R_kF / kF)
            if lamu > 0
                # this condition does not necessarily mean something wrong
                # it only indicates lamu is not converge to correct sign within Ntherm steps
                # normally α>0.8 guarantees convergence, then it means Ntherm is too small
                @warn ("α = $α or Ntherm=$Ntherm is too small!")
            end
            # err = abs(lamu - lamu0)
            # Exit the loop if the iteration converges
            isapprox(lamu, lamu0, rtol=rtol, atol=atol) && break
            n > Nmax && break

            lamu0 = lamu
            if verbose
                println("lamu=$lamu")
            end
        end
    end
    println("α = $α, iteration step: $n")
    lamu = -kF / R0   # calculate 1/R0
    # calculate real physical quantites F and R
    for ni in eachindex(F_freq.mesh[1])
        F_freq[ni, :] ./= kgrid.grid
        R_imt[ni, :] ./= kgrid.grid
    end
    R_ins[1, :] ./= kgrid.grid
    F_fs[:, 1] ./ kF

    println("R(∞)=$(R_ins[ikF])")
    return lamu, F_fs, F_freq, R_imt, R_ins
end

function BSeq_solver_resumB(param, G2::GreenFunc.MeshArray, Πs::GreenFunc.MeshArray,
    kernel, kernel_ins,
    qgrids::Vector{CompositeGrid.Composite},
    Euv; Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8,
    verbose=false, Ncheck=5, Nmax=10000, W=W)

    if verbose
        println("atol=$atol,rtol=$rtol")
    end

    @unpack dim, kF = param
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)
    ikF = kF_label
    iqFs = [locate(qgrids[ki], kF) for ki in 1:kgrid.size]
    println("kF=$kF")
    # println("ikF=$ikF")
    # println("iqFs=$iqFs")

    # B and W in (ω1, ω2, k2)
    B, W0 = initB_resum(W, Euv, G2.mesh[1].representation.rtol, param)
    println("B:$(B.data[1, :, ikF])")
    println("W0:$(W0.data[1, :, ikF])")

    # notice that ω1 is decoupled, and for each ω1 solving B 
    # is exactly iteratively solving R with source W0

    wgrid = B.mesh[1]
    source_ins = GreenFunc.MeshArray([1,], kgrid; dtype=Float64)
    source = GreenFunc.MeshArray(W0.mesh[2], kgrid; dtype=Float64)
    source_ins.data[1, :] .= [-kernel_ins[ik, iqFs[ik]] for ik in 1:length(kgrid)]
    println("source_ins=$(source_ins.data[1,:])")
    for (iw, w) in enumerate(wgrid)
        source.data .= -view(W0.data, iw, :, :)
        println("source(kF)=$(source.data[:,ikF])")
        println("source(0)=$(source.data[1,:])")
        lamu, F_fs, F_freq, R_imt, R_ins = BSeq_solver_resum(param, G2, Πs,
            kernel, kernel_ins,
            qgrids, Euv;
            Ntherm=Ntherm, rtol=rtol, atol=atol, α=α,
            source=source,
            source_ins=source_ins,
            verbose=verbose, Ncheck=Ncheck, Nmax=Nmax)
        Rt2Rw!(F_freq, R_imt, R_ins)
        B.data[iw, :, :] .= F_freq.data
    end
    return B
end

"""
    function linearResponse(param, channel::Int; Euv=100 * param.EF, rtol=1e-10,
        maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8, sigmatype=:none, int_type=:rpa, α=0.7)

Implmentation of Cooper-pair linear response approach.

# Arguments:
- `param`: parameters of ElectronGas.
- `channel::Int`: orbital angular channel (0: s-wave, 1: p-wave, ...)
- `Euv`: the UV energy scale of the spectral density. By default, `Euv=100EF`.
- `rtol`: tolerance absolute error. By defalut, `rtol=1e-10`.
- `maxK`: maximum momentum of kgrid and qgrids. By default, `maxK=10kF`.
- `minK`: minimum interval of panel kgrid and qgrids. By default, `minK=1e-7kF`.
- `Nk`: number of grid points of panel kgrid and qgrids. By defalut, `Nk=8`.
- `order`: number of grid points of subgrid of kgrid and qgrids. By defalut, `order=8`.
- `sigmatype`: type of fermionic self-energy. (no self-energy :none, G0W0 approximation :g0w0)
- `int_type`: type of effective interaction. By default, `int_type=:rpa`.
- `α`: mixing parameter in the self-consistent iteration. By default, `α=0.7`.

# Return:
- Inverse of low-energy linear response `1/R₀` (``R_0=R(\\omega_0, k_F)``)
- Linear response `R(ωₙ, k)`, which is calculated by the Bethe-Slapter-type equation
```math
    R(\\omega_n, k) = 1 - \\frac{1}{\\beta} \\sum_m \\int \\frac{d^dp}{(2\\pi)^d} \\Gamma(\\omega_n,k;\\omega_m,p) G^{(2)}(\\omega_m,p)R(\\omega_m,p)
```
where ``1`` is the default sourced term, ``\\Gamma(k,\\omega_n;p,\\omega_m)`` is the particle-particle four-point vertex 
with zero incoming momentum and frequency, and ``G^{(2)}(p,\\omega_m)`` is the product of two single-particle Green's function.
"""
function linearResponse(param, channel::Int;
    Euv=100 * param.EF, rtol=1e-10, atol=1e-10,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8,
    Vph::Union{Function,Nothing}=nothing, sigmatype=:none, int_type=:rpa,
    α=0.8, verbose=false, Ntherm=30, Nmax=10000,
    resum=false,
    issave=false, uid=1, dir="./", kwargs...)
    println("resum=$resum")
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
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        else
            @time W = LegendreInteraction.DCKernel_old(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
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
        G2 = BSeq.G02wrapped(Euv, rtol, kgrid, param)
    elseif sigmatype == :g0w0
        @unpack me, β, μ = param
        wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
        Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        # self energy should be converted to proper frequency grid
        Σ_dlr = Σ |> to_dlr
        Σ_wn = dlr_to_imfreq(Σ_dlr, wn_mesh)
        G2 = BSeq.G2wrapped(Σ_wn, Σ_ins, param)
    end

    # calculate F, R by Bethe-Slapter iteration.
    if resum
        Πs = Πs0wrapped(Euv, rtol, param)
        lamu, F_fs, F_freq, R_imt, R_ins = BSeq_solver_resum(param, G2, Πs, kernel, kernel_ins, qgrids, Euv;
            rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
        # B = BSeq_solver_resumB(param, G2, Πs, kernel, kernel_ins, qgrids, Euv;
        #     rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax, W=W)
    else
        lamu, F_freq, R_imt, R_ins = BSeq.BSeq_solver(param, G2, kernel, kernel_ins, qgrids, Euv;
            rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
    end
    println("1/R₀ = $lamu")

    R_freq = R_imt |> to_dlr |> to_imfreq
    if issave
        fname = "PCFdata_$(uid).jld2"
        jldopen(dir * fname, "w") do file
            file["param"] = param
            file["lamu"] = lamu
            file["F_freq"] = F_freq
            file["R_ins"] = R_ins
            file["R_freq"] = R_freq
        end
    end
    # println(view(R_freq, :, kF_label))
    return lamu, R_freq, F_freq, R_ins
end


function pcf_resum(param, channel::Int;
    Euv=100 * param.EF, rtol=1e-10, atol=1e-10,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8,
    Vph::Union{Function,Nothing}=nothing, sigmatype=:none, int_type=:rpa,
    α=0.8, verbose=false, Ntherm=30, Nmax=10000,
    issave=false, uid=1, dir="./", kwargs...)
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
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        else
            @time W = LegendreInteraction.DCKernel_old(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        end
        # elseif dim == 2
        #     @time W = LegendreInteraction.DCKernel_2d(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
        #         minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph)
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
        G2 = BSeq.G02wrapped(Euv, rtol, kgrid, param)
    elseif sigmatype == :g0w0
        @unpack me, β, μ = param
        wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
        Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        # self energy should be converted to proper frequency grid
        Σ_dlr = Σ |> to_dlr
        Σ_wn = dlr_to_imfreq(Σ_dlr, wn_mesh)
        G2 = BSeq.G2wrapped(Σ_wn, Σ_ins, param)
    end

    # calculate F, R by Bethe-Slapter iteration.
    Πs = Πs0wrapped(Euv, rtol, param)
    Bwwk = BSeq_solver_resumB(param, G2, Πs, kernel, kernel_ins, qgrids, Euv;
        rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax, W=W)
    lamu, F_fs, F_freq, R_imt, R_ins = BSeq_solver_resum(param, G2, Πs, kernel, kernel_ins, qgrids, Euv;
        rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
    R_freq = R_imt |> to_dlr |> to_imfreq

    A = GreenFunc.MeshArray(R_freq.mesh[1]; dtype=Float64)
    A.data .= R_freq.data[:, kF_label] .+ R_ins.data[1, kF_label]
    B = GreenFunc.MeshArray(Bwwk.mesh[1], Bwwk.mesh[2]; dtype=Float64)
    B.data .= Bwwk.data[:, :, kF_label]

    if issave
        fname = "PCFresum_$(uid).jld2"
        jldopen(dir * fname, "w") do file
            file["param"] = param
            file["A"] = A
            file["B"] = B
        end
    end
    # println(view(R_freq, :, kF_label))
    return A, B
end

end

