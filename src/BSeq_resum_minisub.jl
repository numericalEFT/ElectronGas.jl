# included in BSeq_resum.jl

function calcF_freq!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)
    R_freq = R # assume R_ins is included
    # algebraic in frequency space
    for ind in eachindex(F)
        F[ind] = real(R_freq[ind]) * G2[ind]
    end
end

function initFR_resum_freq_smooth_minisub(wgrid, sgrid, param)
    @unpack β, me, μ = param
    Ω_c = freq_sep

    Rp = GreenFunc.MeshArray(wgrid, sgrid; dtype=Float64)
    Rm = GreenFunc.MeshArray(wgrid, sgrid; dtype=Float64)
    Fp = similar(Rp)
    Fm = similar(Rm)
    for ind in eachindex(Rp)
        ni, ki = ind[1], ind[2]
        ωn, k = wgrid[ni], sgrid[ki]
        e = k^2 / 2 / me - μ
        Rp[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
        Rm[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
    end
    return Fp, Fm, Rp, Rm
end

function initBW_resum_freq_smooth_minisub(W, wgrid, kwgrid, param)
    @unpack dim, kF, β = param

    sgrid = W.kgrid
    kgrid = sgrid
    qgrids = W.qgrids
    bdlr = W.dlrGrid
    ikF = locate(kgrid, kF)
    iqFs = [locate(qgrids[ki], kF) for ki in 1:kgrid.size]

    kernel_ins = W.kernel_bare
    kernel_freq = W.kernel

    @unpack β, me, μ = param
    B = GreenFunc.MeshArray(wgrid, wgrid, [1, 2]; dtype=Float64)

    kernel_dlr = Lehmann.matfreq2dlr(bdlr, kernel_freq; axis=3)
    nec = floor(Int, kwgrid[end] / 2 / π * β + 0.5)
    vn_mesh = GreenFunc.ImFreq(β, BOSON; grid=[i for i in 0:2nec+2])
    # kernel_freq_dense = real(Lehmann.dlr2matfreq(bdlr, kernel_dlr, kwgrid; axis=3))
    kernel_freq_dense = zeros(Float64, (size(kernel_ins, 1), size(kernel_ins, 2), length(kwgrid)))
    for (ik, k) in enumerate(kgrid)
        for (iq, q) in enumerate(qgrids[ik])
            for (iw, w) in enumerate(kwgrid)
                n0 = w / 2 / π * β
                n = Base.floor(Int, n0)
                data = real(Lehmann.dlr2matfreq(bdlr, view(kernel_dlr, ik, iq, :), [n, n + 1]))
                dx0, dx1 = n0 - n, n + 1 - n0
                d0, d1 = data[1], data[2]
                g = d0 * dx1 + d1 * dx0
                gx = g / (dx0 + dx1)
                kernel_freq_dense[ik, iq, iw] = gx
                # if ik == ikF && iq == iqFs[ikF]
                #     println((w, n0, n, data, gx))
                # end
            end
            kernel_freq_dense[ik, iq, :] .+= kernel_ins[ik, iq]
        end
    end
    # println(kernel_freq[ikF, iqFs[ikF], :])
    # kfd = view(kernel_freq_dense, ikF, iqFs[ikF], :) .- kernel_ins[ikF, iqFs[ikF]]
    # println(kfd)
    # println(CompositeGrids.Interp.interp1DGrid(kfd, kwgrid, bdlr.ωn))

    B.data .= 0.0

    return B, kernel_freq_dense
end

function calcR_freq_resum_smooth_minisub!(
    Fp::GreenFunc.MeshArray, Fm::GreenFunc.MeshArray,
    Rp::GreenFunc.MeshArray, Rm::GreenFunc.MeshArray,
    sourcep::GreenFunc.MeshArray, sourcem::GreenFunc.MeshArray,
    ΠB00::Float64,
    kernel, kwgrid, qgrids::Vector{CompositeGrid.Composite})

    kgrid = Fp.mesh[2]
    wgrid = Fp.mesh[1]
    Threads.@threads for ind in eachindex(Rp)
        # for each τ, k, integrate over q
        ωi, ki = ind[1], ind[2]
        # interpolate F to q grid of given k
        qintegrandp = zeros(Float64, length(qgrids[ki]))
        qintegrandm = zeros(Float64, length(qgrids[ki]))
        for (qi, q) in enumerate(qgrids[ki])
            wintegrandp = zeros(Float64, length(wgrid))
            wintegrandm = zeros(Float64, length(wgrid))
            for (wi, w) in enumerate(wgrid)
                wp, wm = wgrid[ωi] + w, abs(wgrid[ωi] - w)
                Fqwp = CompositeGrids.Interp.interp1D(view(Fp, wi, :), kgrid, q)
                Fqwm = CompositeGrids.Interp.interp1D(view(Fm, wi, :), kgrid, q)
                intp = CompositeGrids.Interp.interp1D(view(kernel, ki, qi, :), kwgrid, wp)
                intm = CompositeGrids.Interp.interp1D(view(kernel, ki, qi, :), kwgrid, wm)
                wintegrandp[wi] += (intm * Fqwp + intp * Fqwm)
                wintegrandm[wi] += (intp * Fqwp + intm * Fqwm)
            end
            qintegrandp[qi] += CompositeGrids.Interp.integrate1D(wintegrandp, wgrid)
            qintegrandm[qi] += CompositeGrids.Interp.integrate1D(wintegrandm, wgrid)
        end
        Rp[ωi, ki] = CompositeGrids.Interp.integrate1D(qintegrandp, qgrids[ki]) / 2π / (-4 * π * π)
        Rm[ωi, ki] = CompositeGrids.Interp.integrate1D(qintegrandm, qgrids[ki]) / 2π / (-4 * π * π)
    end
    Rp.data .+= sourcep.data * (1 - ΠB00)
    Rm.data .+= sourcem.data * (1 - ΠB00)
end

function BSeq_solver_freq_resum_smooth_minisub(param, G2::GreenFunc.MeshArray,
    kernel, kwgrid, qgrids::Vector{CompositeGrid.Composite};
    Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8,
    sourcep::GreenFunc.MeshArray=GreenFunc.MeshArray(G2.mesh[1], G2.mesh[2];
        dtype=Float64, data=ones(size(G2))),
    sourcem::GreenFunc.MeshArray=GreenFunc.MeshArray(G2.mesh[1], G2.mesh[2];
        dtype=Float64, data=ones(size(G2))),
    verbose=false, Ncheck=5, Nmax=10000)

    if verbose
        println("atol=$atol,rtol=$rtol")
    end

    @unpack dim, kF = param
    wgrid = G2.mesh[1]
    # iw0 = locate(wgrid.grid, 0)
    iw0 = 1
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)
    ikF = kF_label
    iqFs = [locate(qgrids[ki], kF) for ki in 1:kgrid.size]
    println("kF=$kF")
    # println("ikF=$ikF")
    # println("iqFs=$iqFs")

    for ni in eachindex(sourcep.mesh[1])
        sourcep[ni, :] .*= kgrid.grid
        sourcem[ni, :] .*= kgrid.grid
    end
    source0 = sourcep[iw0, kF_label]
    println("source0=$source0")

    # Initalize F and R
    Fp, Fm, Rp, Rm = initFR_resum_freq_smooth_minisub(wgrid, kgrid, param)

    Rp.data .= sourcep.data
    Rm.data .= sourcem.data

    # R0, R0_sum = 1.0, 0.0
    R0, R0_sum = source0, source0 ./ (1 - α)
    Rp_sum = sourcep.data ./ (1 - α)
    Rm_sum = sourcem.data ./ (1 - α)

    lamu, lamu0 = 0.0, 1.0
    R_kF = Rp[iw0, ikF]
    n = 0
    # self-consistent iteration with mixing α
    # Note!: all quantites about R, F, and source in the `while` loop are multiplied by momentum k.
    while (true)
        n = n + 1
        # switch between imtime and imfreq to avoid convolution
        # dlr Fourier transform is much faster than brutal force convolution

        # calculation from imtime R to imfreq F
        calcF_freq!(Fp, Rp, G2)
        calcF_freq!(Fm, Rm, G2)

        # calculation from imfreq F to imtime R
        if dim == 3 && n > Ntherm
            ω_c = 0.1param.EF
            # ΠB00 = R_kF * (1.0 / 2 / π^2 * log(ω_c * param.β))
            factor = (1 / kF / 4 / π^2 * log(ω_c * param.β / 0.882))
            ΠB00 = Rp[iw0, ikF] * factor
            # ΠB00 = 0.0
            println("ΠB00=$ΠB00, R0=$R0, factor=$(factor)")
            # if abs(ΠB00) > 1.0
            #     ΠB00 = 1.0 * sign(ΠB00)
            # end
            calcR_freq_resum_smooth_minisub!(Fp, Fm, Rp, Rm, sourcep, sourcem, ΠB00, kernel, kwgrid, qgrids)
        end
        R_kF = Rp[iw0, ikF]
        R0_sum = R_kF + R0_sum * α
        R0 = R0_sum * (1 - α)

        # iterative calculation of the dynamical part 
        Rp_sum = view(Rp, :, :) + Rp_sum .* α
        Rp[:, :] = Rp_sum .* (1 - α)
        Rm_sum = view(Rm, :, :) + Rm_sum .* α
        Rm[:, :] = Rm_sum .* (1 - α)
        # @debug "R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))"
        # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

        # record lamu=1/R0 if iterative step n > Ntherm
        if n > Ntherm && (n % Ncheck == 1)
            lamu = -1 / (R_kF / kF)
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
                println("\t R0=$(Rp[iw0,ikF]/kF), R(∞)=$(Rp[end, ikF]/kF)")
                # println("\t R0=$(R0/kF), R(∞)=$(Rp[end, ikF]/kF)")
            end
        end
    end
    println("α = $α, iteration step: $n")
    lamu = -kF / R0   # calculate 1/R0
    # calculate real physical quantites F and R
    for ni in eachindex(Fp.mesh[1])
        Fp[ni, :] ./= kgrid.grid
        Rp[ni, :] ./= kgrid.grid
        Fm[ni, :] ./= kgrid.grid
        Rm[ni, :] ./= kgrid.grid
    end

    println("R0=$(Rp[iw0, ikF])")

    return lamu, Fp, Fm, Rp, Rm
end

function BSeq_solver_resumB_smooth_minisub(param,
    G2::GreenFunc.MeshArray,
    kernel, B, kwgrid,
    qgrids::Vector{CompositeGrid.Composite};
    winit=0, wend=0,
    Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8,
    verbose=false, Ncheck=5, Nmax=10000, W=W,
    issave=false, fname="run/a.jld2")

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

    # notice that ω1 is decoupled, and for each ω1 solving B 
    # is exactly iteratively solving R with source W0

    wgrid = B.mesh[1]
    sourcep = GreenFunc.MeshArray(wgrid, kgrid; dtype=Float64)
    sourcem = GreenFunc.MeshArray(wgrid, kgrid; dtype=Float64)
    if winit == 0
        winit = 1
    end
    if wend == 0
        wend = length(wgrid)
    end
    for iw in winit:wend
        w = wgrid[iw]
        for (ik, k) in enumerate(kgrid)
            for (iv, v) in enumerate(wgrid)
                wp, wm = w + v, abs(w - v)
                intp = CompositeGrids.Interp.interp1D(view(kernel, ik, iqFs[ik], :), kwgrid, wp)
                intm = CompositeGrids.Interp.interp1D(view(kernel, ik, iqFs[ik], :), kwgrid, wm)
                sourcep.data[iv, ik] = -(intm) / k
                sourcem.data[iv, ik] = -(intp) / k
            end
        end
        # println("source(kF)=$(source.data[:,ikF])")
        lamu, Fp, Fm, Rp, Rm = BSeq_solver_freq_resum_smooth_minisub(param, G2,
            kernel, kwgrid, qgrids;
            Ntherm=Ntherm, rtol=rtol, atol=atol, α=α,
            sourcep=sourcep, sourcem=sourcem,
            verbose=verbose, Ncheck=Ncheck, Nmax=Nmax)
        B.data[iw, :, 1] .= Rp.data[:, ikF] ./ kF
        B.data[iw, :, 1] .= Rm.data[:, ikF] ./ kF
        if issave
            jldopen(fname, "w") do file
                file["B"] = B
            end
        end
    end
    return B
end

function pcf_resum_smooth_minisub(param, channel::Int;
    Euv=100 * param.EF, rtol=1e-10, atol=1e-10, Ec=10 * param.EF,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8,
    Vph::Union{Function,Nothing}=nothing, sigmatype=:none, int_type=:rpa,
    α=0.8, verbose=false, Ntherm=30, Nmax=10000,
    onlyA=false, onlyB=false, nB=0, ω_c_ratio=0.1,
    issave=false, uid=1, dir="./", kwargs...)
    @unpack dim, rs, β, kF = param
    if verbose
        println("atol=$atol,rtol=$rtol")
    end
    # prepare Legendre decomposed effective interaction
    if dim == 3
        if channel == 0
            @time W = LegendreInteraction.DCKernel0(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        else
            @time W = LegendreInteraction.DCKernel_old(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        end
    else
        error("No support for $dim dimension!")
    end

    fdlr = Lehmann.DLRGrid(Euv, β, rtol, true, :pha)
    # fdlr = Lehmann.DLRGrid(Euv, β, rtol, true, :none)
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
    fname = "PCFresumdlr_$(uid).jld2"
    alpha = 0.882
    minterval = alpha / β / 2.0
    wgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [alpha / β, Euv], [alpha / β, ω_c_ratio * param.EF], Nk, minterval, order)
    # wgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [-Euv, Euv], [0.0,], Nk, minterval, order)
    kwgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [0.0, 2Euv], [0.0, param.ωp], 16, 0.5minterval, 8)
    # kwgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [0.0, 2Euv], [0.0, param.ωp], 2Nk, minterval, 2order)
    G2 = G02wrapped_freq_smooth(wgrid, kgrid, param)
    # Πs = Πs0wrapped(Euv, rtol, param; ω_c=ω_c_ratio * param.EF)
    # Π = Π0wrapped_freq_smooth(wgrid, param; ω_c=ω_c_ratio * param.EF)
    println(size(G2))
    B, kernel_freq_dense = initBW_resum_freq_smooth_minisub(W, wgrid, kwgrid, param)
    if !(onlyA)
        if nB == 0
            B = BSeq_solver_resumB_smooth_minisub(param, G2, kernel_freq_dense, B, kwgrid, qgrids;
                rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax, W=W, fname=dir * fname, issave=issave)
        else
            B = BSeq_solver_resumB_smooth_minisub(param, G2, kernel_freq_dense, B, kwgrid, qgrids;
                rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax, W=W, fname=dir * fname, issave=issave,
                winit=nB, wend=nB)
        end
    end
    if !(onlyB)
        lamu, F_fs, F_freq, R_freq = BSeq_solver_freq_resum_smooth(param, G2, Πs, kernel_freq_dense, kwgrid, qgrids;
            rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
        A = GreenFunc.MeshArray(R_freq.mesh[1]; dtype=Float64)
        A.data .= R_freq.data[:, kF_label]
    end
    # R_freq = R_imt |> to_dlr |> to_imfreq

    if issave
        jldopen(dir * fname, "w") do file
            file["param"] = param
            if !onlyB
                file["A"] = A
            end
            if !onlyA
                file["B"] = B
            end
        end
    end
    # println(view(R_freq, :, kF_label))
    if onlyA
        return A
    elseif onlyB
        return B
    else
        return A, B
    end

end
