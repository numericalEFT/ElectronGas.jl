
# this file is included in BSeq.jl

# this file provides conventional approach
# solving eigenvalue problem for linearized-Eliashberg equation
# and obtain gap-function at Tc

dot(a, b) = sum(real(a[i] * b[i]) for i in 1:length(a))

function measure_chi(F_freq::GreenFunc.MeshArray)
    F_dlr = F_freq |> to_dlr
    F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
    kgrid = F_ins.mesh[2]
    integrand = view(F_ins, 1, :)
    return real(CompositeGrids.Interp.integrate1D(integrand, kgrid))
end

function eigen_solver(param,
    G2::GreenFunc.MeshArray, kernel, kernel_ins,
    qgrids::Vector{CompositeGrid.Composite}, Euv;
    shift=2.0,
    Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8,
    verbose=false, Ncheck=5, Nmax=10000)

    if verbose
        println("atol=$atol,rtol=$rtol")
    end

    @unpack dim, kF = param
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)

    # Initalize F and R
    F_freq, R_imt, R_ins = initFR(Euv, G2.mesh[1].representation.rtol, kgrid, param)
    R0old = similar(R_ins.data)
    Rold = similar(R_imt.data)
    # R0old .= 1.0
    # Rold .= 0.0

    lamu, lamu0 = 0.0, 1.0
    n = 0
    # self-consistent iteration with mixing α
    # Note!: all quantites about R, F, and source in the `while` loop are multiplied by momentum k.
    while (true)
        n = n + 1
        R0old .= R_ins.data
        Rold .= R_imt.data

        # switch between imtime and imfreq to avoid convolution
        # dlr Fourier transform is much faster than brutal force convolution

        # calculation from imtime R to imfreq F
        calcF!(F_freq, R_imt, R_ins, G2)

        # calculation from imfreq F to imtime R
        if dim == 3
            calcR!(F_freq, R_imt, R_ins, Nothing(), kernel, kernel_ins, qgrids)
        elseif dim == 2
            calcR_2d!(F_freq, R_imt, R_ins, Nothing(), kernel, kernel_ins, qgrids)
        else
            error("Not implemented for $dim dimensions.")
        end

        lamu = dot(R_ins.data, R0old)

        R_ins[:, :] = R_ins[:, :] .+ shift .* R0old
        R_imt[:, :] = R_imt[:, :] .+ shift .* Rold

        modulus = sqrt(dot(R_ins.data, R_ins.data))

        R_ins[:, :] = R_ins[:, :] ./ modulus
        R_imt[:, :] = R_imt[:, :] ./ modulus

        # println("$lamu, $modulus")
        # record lamu if iterative step n > Ntherm
        if n > Ntherm && (n % Ncheck == 1)
            isapprox(lamu, lamu0, rtol=rtol, atol=atol) && break
            n > Nmax && break

            lamu0 = lamu
            if verbose
                println("lamu=$lamu")
            end
        end
    end
    println("iteration step: $n")
    # calculate real physical quantites F and R
    for ni in eachindex(F_freq.mesh[1])
        F_freq[ni, :] ./= kgrid.grid
        R_imt[ni, :] ./= kgrid.grid
    end
    R_ins[1, :] ./= kgrid.grid

    return lamu, F_freq, R_imt, R_ins
end

function lin_eliashberg(param, channel::Int;
    shift=2.0,
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
    lamu, F_freq, R_imt, R_ins = eigen_solver(param, G2, kernel, kernel_ins, qgrids, Euv; shift=shift,
        rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
    println("lamu = $lamu")

    R_freq = R_imt |> to_dlr |> to_imfreq
    if issave
        fname = "gapdata_$(uid).jld2"
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
