"""
Calculate gap-function equation
"""

using LinearAlgebra
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids
using ElectronGas.Convention

const dim = 2
# const beta, rs = 1e5, 1.5

function G02wrapped(fdlr, kgrid, param)
    # return G0(K)G0(-K)
    @unpack me, kF, β, μ = param

    green_dyn = zeros(Float64, (kgrid.size, fdlr.size))
    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(fdlr.n)
            ω = k^2 / 2 / me - μ
            green_dyn[ki, ni] = 1 / (
                ((2n + 1) * π / β)^2
                +
                (ω)^2
            )
        end
    end

    return green_dyn
end

function G2wrapped(fdlr, kgrid, param, Σ::GreenFunc.Green2DLR)
    # return G(K)G(-K)
    @unpack me, kF, β, μ = param

    Σ_freq = GreenFunc.toMatFreq(Σ, fdlr.n)
    Σ_shift = real(GreenFunc.dynamic(Σ_freq, π / β, kF, 1, 1) + GreenFunc.instant(Σ_freq, kF, 1, 1))

    green_dyn = zeros(Float64, (kgrid.size, fdlr.size))

    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(fdlr.n)
            ω = k^2 / 2 / me - μ
            ΣR, ΣI = real(Σ_freq.dynamic[1, 1, ki, ni] + Σ_freq.instant[1, 1, ki] - Σ_shift), imag(Σ_freq.dynamic[1, 1, ki, ni])
            green_dyn[ki, ni] = 1 / (
                ((2n + 1) * π / β - ΣI)^2
                +
                (ω + ΣR)^2
            )
        end
    end
    return green_dyn
end

function ΔFinit(fdlr, kgrid)

    delta = zeros(Float64, (kgrid.size, fdlr.size))
    F = zeros(Float64, (kgrid.size, fdlr.size))
    delta0 = zeros(Float64, kgrid.size)

    for (ki, k) in enumerate(kgrid)
        delta0[ki] = 1.0
        for (τi, τ) in enumerate(fdlr.τ)
            delta[ki, τi] = 1.0
        end
    end

    return delta, delta0, F
end

function dotΔ(fdlr, kgrid, Δ, Δ0, Δ2=Δ, Δ02=Δ0)
    # kF_label = searchsortedfirst(kgrid.grid, param.kF)
    # return Δ0[kF_label] + real(Lehmann.tau2matfreq(fdlr, view(Δ, kF_label, :), fdlr.n))[1]
    # return (dot(Δ, Δ2) + dot(Δ0, Δ02))
    return (dot(Δ, Δ2))
end

function calcF!(F, fdlr, kgrid, Δ, Δ0, G2)
    Δ_freq = real(Lehmann.tau2matfreq(fdlr, Δ, fdlr.n; axis=2))
    for (ki, k) in enumerate(kgrid.grid)
        for (ni, n) in enumerate(fdlr.n)
            F[ki, ni] = (Δ_freq[ki, ni] + Δ0[ki]) * G2[ki, ni]

        end
    end
end

function calcΔ!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)

    F_tau = real(Lehmann.matfreq2tau(fdlr, F, fdlr.τ; axis=2))
    F_ins = -real(tau2tau(fdlr, F_tau, [fdlr.β,]; axis=2))[:, 1]
    for (ki, k) in enumerate(kgrid.grid)
        for (τi, τ) in enumerate(fdlr.τ)
            Fk = CompositeGrids.Interp.interp1DGrid(view(F_tau, :, τi), kgrid, qgrids[ki].grid)
            integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fk
            Δ[ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            if τi == 1
                Fk = CompositeGrids.Interp.interp1DGrid(F_ins, kgrid, qgrids[ki].grid)
                integrand = view(kernel_bare, ki, 1:qgrids[ki].size) .* Fk
                Δ0[ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            end
        end
    end
end

function calcΔ_2d!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)

    F_tau = real(Lehmann.matfreq2tau(fdlr, F, fdlr.τ; axis=2))
    F_ins = -real(tau2tau(fdlr, F_tau, [fdlr.β,]; axis=2))[:, 1]
    for (ki, k) in enumerate(kgrid.grid)
        for (τi, τ) in enumerate(fdlr.τ)
            Fk = CompositeGrids.Interp.interp1DGrid(view(F_tau, :, τi), kgrid, qgrids[ki].grid)
            integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fk .* k
            Δ[ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            if τi == 1
                Fk = CompositeGrids.Interp.interp1DGrid(F_ins, kgrid, qgrids[ki].grid)
                integrand = view(kernel_bare, ki, 1:qgrids[ki].size) .* Fk .* k
                Δ0[ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            end
        end
    end
end

function gapIteration(param, fdlr, kgrid, qgrids, kernel, kernel_bare, G2;
    Nstep=1e2, rtol=1e-6, shift=2.0)
    @unpack dim = param

    Δ, Δ0, F = ΔFinit(fdlr, kgrid)

    delta = zeros(Float64, (kgrid.size, fdlr.size))
    delta0 = zeros(Float64, (kgrid.size))

    n = 0
    lamu, lamu0 = 1.0, 2.0
    err = 1.0

    while (n < Nstep && err > rtol)

        calcF!(F, fdlr, kgrid, Δ, Δ0, G2)

        n = n + 1

        delta = copy(Δ)
        delta0 = copy(Δ0)

        if dim == 3
            calcΔ!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)
        elseif dim == 2
            calcΔ_2d!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)
        end

        lamu = dotΔ(fdlr, kgrid, Δ, Δ0, delta, delta0)
        # dotΔ(fdlr, kgrid, Δ, Δ0, Δ2 = Δ, Δ02 = Δ0)

        Δ0 = Δ0 .+ shift .* delta0
        Δ = Δ .+ shift .* delta

        #modulus = Normalization(delta_0_new, delta_0_new, kgrid)
        modulus = sqrt(dotΔ(fdlr, kgrid, Δ, Δ0))

        Δ = Δ ./ modulus
        Δ0 = Δ0 ./ modulus
        err = abs(lamu - lamu0) / abs(lamu + EPS)
        lamu0 = lamu
        # println(lamu)
    end
    return lamu, Δ, Δ0, F
end



# if abspath(PROGRAM_FILE) == @__FILE__
function gapfunction(beta, rs,
    channel::Int, dim::Int, sigmatype;
    issave=false, uid=100)
    #--- parameters ---
    param = Parameter.defaultUnit(1 / beta, rs, dim)
    # Euv, rtol = 100 * param.EF, 1e-11
    # maxK, minK = 20param.kF, 1e-8param.kF
    # Nk, order = 12, 10
    Euv, rtol = 100 * param.EF, 1e-10
    maxK, minK = 10param.kF, 1e-7param.kF
    Nk, order = 8, 8
    int_type = :rpa
    # print(param)

    #--- prepare kernel ---
    if dim == 3
        @time W = LegendreInteraction.DCKernel0(param;
            Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
            int_type=int_type, channel=channel)
    elseif dim == 2
        @time W = LegendreInteraction.DCKernel_2d(param;
            Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
            int_type=int_type, channel=channel)
    else
        error("No support for $dim dimension!")
    end

    fdlr = Lehmann.DLRGrid(Euv, param.β, rtol, true, :pha)
    bdlr = W.dlrGrid
    kgrid = W.kgrid
    qgrids = W.qgrids

    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3))

    kF_label = searchsortedfirst(kgrid.grid, param.kF)
    qF_label = searchsortedfirst(qgrids[kF_label].grid, param.kF)

    println("static kernel at (kF, kF):$(kernel_bare[kF_label, qF_label])")
    # println("dynamic kernel at (kF, kF):")
    # println(view(kernel_freq, kF_label, qF_label, :))

    #--- prepare Σ ---
    if sigmatype == :g0w0
        Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
    end

    #--- prepare G2 ---
    if sigmatype == :none
        G2 = G02wrapped(fdlr, kgrid, param)
    elseif sigmatype == :g0w0
        G2 = G2wrapped(fdlr, kgrid, param, Σ)
    end

    # println("G2 at kF:")
    # println(view(G2, kF_label, :))

    lamu, Δ, Δ0, F = gapIteration(param, fdlr, kgrid, qgrids, kernel, kernel_bare, G2; shift=4.0)
    println("beta = $beta,    rs = $rs")
    println("channel = $channel")
    println("lamu = $lamu")

    Δ_freq = real(Lehmann.tau2matfreq(fdlr, Δ, fdlr.n; axis=2))
    println(fdlr.ωn ./ param.EF)
    println(view(Δ_freq, kF_label, :))
    # println(view(F, kF_label, :))

    if issave
        fname = "gapdata_$(uid).jld2"
        jldopen(dir * fname, "w") do file
            file["param"] = param
            file["lamu"] = lamu
            file["F_freq"] = F
            file["R_ins"] = Δ0
            file["R_freq"] = Δ_freq
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    sigmatype = :none
    # sigmatype = :g0w0

    rs = 2.0
    # rs = 0.5
    channel = 0
    # blist = [1e2, 1e3, 1e4, 1e5, 2e5]
    blist = [1e2]
    # blist = [1e2, 1e3, 5e3, 1e4]

    for (ind, beta) in enumerate(blist)
        gapfunction(beta, rs, channel, dim, sigmatype)
    end
end