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

function G02wrapped(fdlr, kgrid, param)
    # return G(K)G(-K)
    @unpack me, kF, β, EF = param

    green_dyn = zeros(Float64, (kgrid.size, fdlr.size))
    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(fdlr.n)
            ω = k^2/2/me
            green_dyn[ki,ni] = 1/(
                ( (2n + 1) * π / β) ^ 2
                + (ω) ^ 2
            )
        end
    end

    return green_dyn
end

function ΔFinit(fdlr, kgrid)
    @unpack me, kF, β, EF = param

    delta = zeros(Float64, (kgrid.size, fdlr.size))
    F = zeros(Float64, (kgrid.size, fdlr.size))
    delta0 = zeros(Float64, kgrid.size)

    for (ki, k) in enumerate(kgrid)
        delta0[ki] = 1.0
        for (τi, τ) in enumerate(fdlr.τ)
            delta[ki, τi] = 0.0
        end
    end

    return delta, delta0, F
end

function dotΔ(fdlr, kgrid, Δ, Δ0, Δ2 = Δ, Δ02 = Δ0)
    # kF_label = searchsortedfirst(kgrid.grid, param.kF)
    # return Δ0[kF_label] + real(Lehmann.tau2matfreq(fdlr, view(Δ, kF_label, :), fdlr.n))[1]
    return (dot(Δ, Δ2) + dot(Δ0, Δ02))
end

function calcF!(F, fdlr, kgrid, Δ, Δ0, G2)
    Δ_freq = real(Lehmann.tau2matfreq(fdlr, Δ, fdlr.n; axis = 2))
    for (ki, k) in enumerate(kgrid.grid)
        for (ni, n) in enumerate(fdlr.n)
            F[ki, ni] = (Δ_freq[ki, ni] + Δ0[ki]) * G2[ki, ni]

        end
    end
end

function calcΔ!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)

    F_tau = real(Lehmann.matfreq2tau(fdlr, F, fdlr.τ; axis = 2))
    F_ins = - real(tau2tau(fdlr, F_tau, [fdlr.β,]; axis = 2))[:, 1]
    for (ki, k) in enumerate(kgrid.grid)
        for (τi, τ) in enumerate(fdlr.τ)
            Fk = CompositeGrids.Interp.interp1DGrid(view(F_tau, :, τi), kgrid, qgrids[ki].grid)
            integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fk
            Δ[ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])./(-4*π*π)
            if τi == 1
                Fk = CompositeGrids.Interp.interp1DGrid(F_ins, kgrid, qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Fk
                Δ0[ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])./(-4*π*π)
            end
        end
    end
end

function gapIteration(param, fdlr, kgrid, qgrids,  kernel, kernel_bare, G2;
                      Nstep = 1e2, rtol = 1e-5, shift = 2.0)

    Δ, Δ0, F = ΔFinit(fdlr, kgrid)

    delta = zeros(Float64, (kgrid.size, fdlr.size))
    delta0 = zeros(Float64, (kgrid.size))

    n = 0
    lamu, lamu0 = 1.0, 2.0
    err = 1.0

    while(n < Nstep && err > rtol)

        calcF!(F, fdlr, kgrid, Δ, Δ0, G2)

        n=n+1

        delta = copy(Δ)
        delta0 = copy(Δ0)

        calcΔ!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)

        lamu = dotΔ(fdlr, kgrid, Δ, Δ0, delta, delta0) ^ 2

        Δ0 = Δ0 .+ shift .* delta0
        Δ = Δ .+ shift .* delta

        #modulus = Normalization(delta_0_new, delta_0_new, kgrid)
        modulus = sqrt(dotΔ(fdlr, kgrid, Δ, Δ0))

        Δ = Δ ./ modulus
        Δ0 = Δ0 ./ modulus
        err=abs(lamu-lamu0)/abs(lamu+EPS)
        lamu0=lamu
        println(lamu)
    end
    return lamu, Δ, Δ0, F
end



if abspath(PROGRAM_FILE) == @__FILE__

    #--- parameters ---


    param = Parameter.defaultUnit(1/250.0, 4.0)
    Euv, rtol = 100*param.EF, 1e-10
    maxK, minK = 10param.kF, 1e-7param.kF
    Nk, order = 16, 4
    int_type = :rpa

    #--- prepare kernel ---
    W = LegendreInteraction.DCKernel0(param;
                                      Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order,
                                      int_type = int_type)

    fdlr = Lehmann.DLRGrid(Euv, param.β, rtol, true, :pha)
    bdlr = W.dlrGrid
    kgrid = W.kgrid
    qgrids = W.qgrids

    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis = 3))

    kF_label = searchsortedfirst(kgrid.grid, param.kF)
    qF_label = searchsortedfirst(qgrids[kF_label].grid, param.kF)

    println("dynamic kernel at (kF, kF):")
    println(view(kernel, kF_label, qF_label, :))

    #--- prepare G2 ---

    G2 = G02wrapped(fdlr, kgrid, param)

    println("G2 at kF:")
    println(view(G2, kF_label, :))

    #--- prepare Δ ---
    Δ, Δ0, F = ΔFinit(fdlr, kgrid)
    norm = dotΔ(fdlr, kgrid, Δ, Δ0)
    println("norm=$norm")
    println(view(Δ, kF_label, :))

    calcF!(F, fdlr, kgrid, Δ, Δ0, G2)
    println(view(F, kF_label, :))

    calcΔ!(Δ, Δ0, fdlr, kgrid, qgrids, F, kernel, kernel_bare)
    println(view(Δ, kF_label, :))

    lamu, Δ, Δ0, F = gapIteration(param, fdlr, kgrid, qgrids, kernel, kernel_bare, G2)
    println("lamu=$lamu")
    println(view(Δ, kF_label, :))
    println(view(F, kF_label, :))

end
