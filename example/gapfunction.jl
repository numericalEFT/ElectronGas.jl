"""
Calculate gap-function equation
"""

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids

function G2wrapped(param, Σ::GreenFunc.Green2DLR)
    # return G(K)G(-K)
    @unpack me, kF, β, EF = param
    Σ_freq = GreenFunc.toMatFreq(Σ)
    Σ_shift = real(GreenFunc.dynamic(Σ_freq, π / β, kF, 1, 1) + GreenFunc.instant(Σ_freq, kF, 1, 1))
    green = Green2DLR{Float64}(
        :G, GreenFunc.IMFREQ, Σ_freq.β, Σ_freq.isFermi, Σ_freq.dlrGrid.Euv, Σ_freq.spaceGrid, Σ_freq.color;
        timeSymmetry = Σ_freq.timeSymmetry, rtol = Σ_freq.dlrGrid.rtol)

    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(green.spaceGrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            ω = k^2 / 2 / me
            ΣR, ΣI = real(Σ_freq.dynamic[1, 1, ki, ni] + Σ_freq.instant[1, 1, ki] - Σ_shift), imag(Σ_freq.dynamic[1, 1, ki, ni])
            green_dyn[1, 1, ki, ni] = 1 / (
                ((2n + 1) * π / β - ΣI)^2
                +
                (ω + ΣR)^2
            )
        end
    end
    green.dynamic = green_dyn
    return green
end

function G02wrapped(param, Euv, rtol, kgrid)
    # return G(K)G(-K)
    @unpack me, kF, β, EF = param
    green = Green2DLR{Float64}(
        :G, GreenFunc.IMFREQ, β, true, Euv, kgrid, 1;
        timeSymmetry = :pha, rtol = rtol)

    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(green.spaceGrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            ω = k^2 / 2 / me
            green_dyn[1, 1, ki, ni] = 1 / (
                ((2n + 1) * π / β)^2
                +
                (ω)^2
            )
        end
    end
    green.dynamic = green_dyn
    return green
end

function normalize!(Δ::GreenFunc.Green2DLR)
    kgrid = Δ.spaceGrid
    λ = real(CompositeGrids.Interp.integrate1D(Δ.instant[1, 1, :] .+ Δ.dynamic[1, 1, :, 1], kgrid))
    Δ.instant = Δ.instant ./ λ
    Δ.dynamic = Δ.dynamic ./ λ
    return λ
end

function calcΔ!(F::GreenFunc.Green2DLR, W::LegendreInteraction.DCKernel, Δ::GreenFunc.Green2DLR)
    @unpack β = W.param
    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = F.dlrGrid
    bdlr = W.dlrGrid

    kernel_bare = W.kernel_bare ./ (4 * π * π)
    kernel_freq = W.kernel
    kernel = Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis = 3) ./ (4 * π * π)

    F_tau = GreenFunc.toTau(F)
    F_ins = real(tau2tau(F_tau.dlrGrid, F_tau.dynamic, [β,], F_tau.timeGrid.grid; axis = 4)[1, 1, :, 1])

    for (ki, k) in enumerate(Δ.spaceGrid)
        for (τi, τ) in enumerate(Δ.dlrGrid.τ)
            Fk = CompositeGrids.Interp.interp1DGrid(F.dynamic[1, 1, :, τi], kgrid, qgrids[ki].grid)
            integrand = kernel[ki, 1:qgrids[ki].size, τi] .* Fk ./ k .* qgrids[ki].grid
            Δ.dynamic[1, 1, ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            if τi == 1
                Fk = CompositeGrids.Interp.interp1DGrid(F_ins, kgrid, qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Fk ./ k .* qgrids[ki].grid
                Δ.instant[1, 1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            end
        end
    end
end

function calcF!(Δ::GreenFunc.Green2DLR, G2::GreenFunc.Green2DLR, F::GreenFunc.Green2DLR)
    Δ_freq = GreenFunc.toMatFreq(Δ)
    for (ki, k) in enumerate(F.spaceGrid)
        for (ni, n) in enumerate(F.dlrGrid.n)
            F.dynamic[1, 1, ki, ni] = real(Δ_freq.dynamic[1, 1, ki, ni] + Δ_freq.instant[1, 1, ki]) * G2.dynamic[1, 1, ki, ni]
        end
    end
end

function gapIteration(param, N, Euv, rtol, Nk, maxK, minK, order, int_type)
    @unpack β = param

    W = LegendreInteraction.DCKernel0(param;
        Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type)
    kgrid = W.kgrid
    qgrids = W.qgrids
    bdlr = W.dlrGrid

    # Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
    G2 = G02wrapped(param, Euv, rtol, W.kgrid)
    fdlr = G2.dlrGrid


    Δ = GreenFunc.Green2DLR{Float64}(
        :delta, GreenFunc.IMTIME, param.β, true, G2.dlrGrid.Euv, kgrid, 1;
        timeSymmetry = :pha, rtol = rtol)
    Δ.instant = ones(Float64, (1, 1, length(kgrid.grid)))
    Δ.dynamic = ones(Float64, (1, 1, length(kgrid.grid), fdlr.size))

    F = GreenFunc.Green2DLR{Float64}(
        :F, GreenFunc.IMFREQ, param.β, true, G2.dlrGrid.Euv, W.kgrid, 1;
        timeSymmetry = :pha, rtol = rtol)
    F.instant = ones(Float64, (1, 1, length(kgrid.grid)))
    F.dynamic = ones(Float64, (1, 1, length(kgrid.grid), fdlr.size))

    Δ_freq = GreenFunc.toMatFreq(Δ)

    shift = 2.0

    for i in 1:N
        Δ_freq_new = GreenFunc.toMatFreq(Δ)
        Δ_freq.instant = Δ_freq.instant .* shift .+ Δ_freq_new.instant
        Δ_freq.dynamic = Δ_freq.dynamic .* shift .+ Δ_freq_new.dynamic

        λ = normalize!(Δ_freq) - shift
        println("λ = $λ")

        calcF!(Δ_freq, G2, F)
        calcΔ!(F, W, Δ)
    end

    return Δ
end

if abspath(PROGRAM_FILE) == @__FILE__
    param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(1 / 1000.0, 3.0)
    Euv, rtol = 100 * param.EF, 1e-10
    Nk, order = 8, 4

    Δ = gapIteration(param, 10, Euv, rtol, Nk, 10 * param.kF, 1e-7 * param.kF, order, :rpa)

    kgrid = Δ.spaceGrid
    kF = kgrid.panel[3]
    kF_label = searchsortedfirst(kgrid.grid, kF)

    println(Δ.instant[1, 1, :])
    println(Δ.dynamic[1, 1, kF_label, :])
end
