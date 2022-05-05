"""
Functions in the gap-equation solver.
"""

using LinearAlgebra
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids
using ElectronGas.Convention


function calcF!(F::GreenFunc.Green2DLR, Δ::GreenFunc.Green2DLR, G2::GreenFunc.Green2DLR)
    Δ_freq = tau2matfreq(Δ.dlrGrid, Δ.dynamic, Δ.dlrGrid.n, Δ.timeGrid.grid; axis=4)
    # println(real(Δ_freq.dynamic[1, 1, 113, :]))
    # println(real(Δ_freq.instant[1, 1, 113]))
    # println(real(Δ.instant[1, 1, 113]))

    for (ki, k) in enumerate(F.spaceGrid.grid)
        for (ni, n) in enumerate(F.dlrGrid.n)
            F.dynamic[1, 1, ki, ni] = real(Δ_freq[1, 1, ki, ni] + Δ.instant[1, 1, ki]) * G2.dynamic[1, 1, ki, ni]
        end
    end
end

function calcΔ!(F::GreenFunc.Green2DLR, Δ::GreenFunc.Green2DLR, kernel, kernel_bare, qgrids)
    kgrid = F.spaceGrid

    F_tau = real(matfreq2tau(F.dlrGrid, F.dynamic, F.dlrGrid.τ, F.timeGrid.grid; axis=4))
    F_ins = -real(tau2tau(F.dlrGrid, F_tau, [F.dlrGrid.β,], F.dlrGrid.τ; axis=4))[1, 1, :, 1]
    for (ki, k) in enumerate(kgrid.grid)
        for (τi, τ) in enumerate(Δ.dlrGrid.τ)
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_tau, 1, 1, :, τi), kgrid, qgrids[ki].grid)
            integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq
            Δ.dynamic[1, 1, ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            if τi == 1
                Fq = CompositeGrids.Interp.interp1DGrid(F_ins, kgrid, qgrids[ki].grid)
                integrand = view(kernel_bare, ki, 1:qgrids[ki].size) .* Fq
                # Δ.instant[1, 1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
                Δ.instant[1, 1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            end
        end
    end
end

function calcΔ_2d!(F::GreenFunc.Green2DLR, Δ::GreenFunc.Green2DLR, kernel, kernel_bare, qgrids)
    kgrid = F.spaceGrid

    F_tau = real(matfreq2tau(F.dlrGrid, F.dynamic, F.dlrGrid.τ, F.timeGrid.grid; axis=4))
    F_ins = -real(tau2tau(F.dlrGrid, F_tau, [F.dlrGrid.β,], F.dlrGrid.τ; axis=4))[1, 1, :, 1]
    for (ki, k) in enumerate(kgrid.grid)
        for (τi, τ) in enumerate(Δ.dlrGrid.τ)
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_tau, 1, 1, :, τi), kgrid, qgrids[ki].grid)
            integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq .* k
            Δ.dynamic[1, 1, ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            if τi == 1
                Fq = CompositeGrids.Interp.interp1DGrid(F_ins, kgrid, qgrids[ki].grid)
                integrand = view(kernel_bare, ki, 1:qgrids[ki].size) .* Fq .* k
                Δ.instant[1, 1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            end
        end
    end
end