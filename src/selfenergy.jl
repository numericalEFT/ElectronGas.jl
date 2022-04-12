"""
Calculate self-energy
"""

module SelfEnergy

using ..Parameter, ..Convention, ..Polarization, ..Interaction, ..LegendreInteraction
using ..Parameters, ..GreenFunc, ..Lehmann, ..LegendrePolynomials, ..CompositeGrids

@inline function Fock0_3dZeroTemp(k, param)
    @assert param.ga ≈ 0 "current implementation only supports spin-symmetric interaction"
    @assert param.dim == 3
    @assert k >= 0
    # TODO: add spin-asymmetric interaction
    @unpack me, kF, Λs, e0 = param

    if k < 1e-6
        k = 1e-6
    end

    l = sqrt(Λs)
    if l > 1e-14
        fock = 1 + l / kF * (atan((k - kF) / l) - atan((k + kF) / l))
        fock -= (l^2 - k^2 + kF^2) / 4 / k / kF * log((l^2 + (k - kF)^2) / (l^2 + (k + kF)^2))
    else
        if abs(k - kF) > 1e-10
            fock = 1 + (k^2 - kF^2) / 4 / k / kF * log((k - kF)^2 / (k + kF)^2)
        else
            fock = 1
        end
    end

    return fock * (-e0^2 * kF) / π

end

@inline function Fock0_2dZeroTemp(k, param)
    @assert param.ga ≈ 0 "current implementation only supports spin-symmetric interaction"
    @assert param.dim == 2
    # TODO: add spin-asymmetric interaction
    @unpack me, kF, Λs, e0 = param
    l2 = Λs
    x = √((k - kF)^2 + l2)
    y = √((k + kF)^2 + l2)
    z = √(k^2 + l2)

    if abs(k) > EPS
        if l2 < EPS
            fock = abs(k - kF) * (1 - kF / k) + abs(k + kF) * (1 + kF / k) - 2 * abs(k)
        else
            fock = (k - kF) * x / k + (k + kF) * y / k - 2z
            fock += l2 * log((k - kF + x) * (k + kF + y) / (k + z)^2) / k
        end
    else
        fock = 4 * (√(kF^2 + l2) - sqrt(Λs))
    end

    return fock * (-e0^2) / 4π
end

# @inline function Fock0_2dZeroTemp(k, param)
#     @assert param.e0a ≈ 0 "current implementation only supports spin-symmetric interaction"
#     @assert param.dim == 2
#     # TODO: add spin-asymmetric interaction
#     @unpack me, kF, Λs, e0 = param
#     @assert !(Λs ≈ 0.0) "Fock diverges diverges for the bare Coulomb interaction"
#     l2 = Λs
#     x = kF^2 + l2 - k^2
#     c = 4 * k^2 * l2
#     return -e0^2 * log((sqrt(x^2 + c) + x) / 2 / l2)
# end

"""
    function Fock0_ZeroTemp(q, n, param)

Zero temperature one-spin Fock function for momentum.
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F) and Yukawa/Coulomb instant interaction.


#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
  - param: other system parameters
"""
function Fock0_ZeroTemp(k::Float64, param)
    @unpack dim = param

    if dim == 2
        return Fock0_2dZeroTemp(k, param)
    elseif dim == 3
        return Fock0_3dZeroTemp(k, param)
    else
        error("No support for zero-temperature Fock in $dim dimension!")
    end
end

function G0wrapped(Euv, rtol, sgrid, param)
    @unpack me, kF, β, μ = param

    green = GreenFunc.Green2DLR{ComplexF64}(:g0, GreenFunc.IMFREQ, β, true, Euv, sgrid, 1; rtol=rtol)
    green_dyn = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1, 1, ki, ni] = 1 / (im * (π / β * (2n + 1)) - (k^2 / 2 / me - μ))
        end
    end
    green.dynamic = green_dyn
    return green
end

function Gwrapped(Σ::GreenFunc.Green2DLR, param)
    @unpack me, kF, β, μ = param
    Σ_freq = GreenFunc.toMatFreq(Σ)
    green = Green2DLR{ComplexF64}(
        :G, GreenFunc.IMFREQ, Σ_freq.β, Σ_freq.isFermi, Σ_freq.dlrGrid.Euv, Σ_freq.spaceGrid, Σ_freq.color;
        timeSymmetry=Σ_freq.timeSymmetry, rtol=Σ_freq.dlrGrid.rtol)

    green_dyn = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    Σ_shift = real(GreenFunc.dynamic(Σ_freq, π / β, kF, 1, 1) + GreenFunc.instant(Σ_freq, kF, 1, 1))

    for (ki, k) in enumerate(green.spaceGrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1, 1, ki, ni] = 1 / (im * (π / β * (2n + 1)) - (k^2 / 2 / me - μ) -
                                           Σ.dynamic[1, 1, ki, ni] - Σ.instant[1, 1, ki] + Σ_shift)
        end
    end
    green.dynamic = green_dyn
    return green
end

function calcΣ_2d(G::GreenFunc.Green2DLR, W::LegendreInteraction.DCKernel)
    @unpack β = W.param

    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = G.dlrGrid
    bdlr = W.dlrGrid
    G = GreenFunc.toTau(G)

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3)

    # container of Σ
    Σ = GreenFunc.Green2DLR{ComplexF64}(:sigma, GreenFunc.IMTIME, β, true, fdlr.Euv, kgrid, 1; rtol=fdlr.rtol)
    Σ_ins = zeros(ComplexF64, (1, 1, length(kgrid.grid)))
    Σ_dyn = zeros(ComplexF64, (1, 1, length(kgrid.grid), fdlr.size))

    # equal-time green (instant)
    G_ins = tau2tau(G.dlrGrid, G.dynamic, [β,], G.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
    @assert length(G.instant) == 0 "current implication supports green function without instant part"

    for (ki, k) in enumerate(kgrid.grid)

        for (τi, τ) in enumerate(fdlr.τ)
            Gq = CompositeGrids.Interp.interp1DGrid(G.dynamic[1, 1, :, τi], kgrid, qgrids[ki].grid)
            integrand = kernel[ki, 1:qgrids[ki].size, τi] .* Gq .* qgrids[ki].grid
            Σ_dyn[1, 1, ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            @assert isfinite(Σ_dyn[1, 1, ki, τi]) "fail Δ at $ki, $τi"
            if τi == 1
                Gq = CompositeGrids.Interp.interp1DGrid(G_ins, kgrid, qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Gq .* qgrids[ki].grid
                Σ_ins[1, 1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
                @assert isfinite(Σ_ins[1, 1, ki]) "fail Δ0 at $ki"
            end
        end
    end

    Σ.dynamic, Σ.instant = Σ_dyn ./ (-4 * π^2), Σ_ins ./ (-4 * π^2)
    return Σ
end

# function calcΣ(kernal, kernal_bare, fdlr, kgrid, qgrids)
function calcΣ_3d(G::GreenFunc.Green2DLR, W::LegendreInteraction.DCKernel)
    @unpack β = W.param

    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = G.dlrGrid
    bdlr = W.dlrGrid
    G = GreenFunc.toTau(G)

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3)

    # container of Σ
    Σ = GreenFunc.Green2DLR{ComplexF64}(:sigma, GreenFunc.IMTIME, β, true, fdlr.Euv, kgrid, 1; rtol=fdlr.rtol)
    Σ_ins = zeros(ComplexF64, (1, 1, length(kgrid.grid)))
    Σ_dyn = zeros(ComplexF64, (1, 1, length(kgrid.grid), fdlr.size))

    # equal-time green (instant)
    G_ins = tau2tau(G.dlrGrid, G.dynamic, [β,], G.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
    @assert length(G.instant) == 0 "current implication supports green function without instant part"

    for (ki, k) in enumerate(kgrid.grid)
        for (τi, τ) in enumerate(fdlr.τ)
            Gq = CompositeGrids.Interp.interp1DGrid(G.dynamic[1, 1, :, τi], kgrid, qgrids[ki].grid)
            integrand = kernel[ki, 1:qgrids[ki].size, τi] .* Gq ./ k .* qgrids[ki].grid
            Σ_dyn[1, 1, ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            @assert isfinite(Σ_dyn[1, 1, ki, τi]) "fail Δ at $ki, $τi"
            if τi == 1
                Gq = CompositeGrids.Interp.interp1DGrid(G_ins, kgrid, qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Gq ./ k .* qgrids[ki].grid
                Σ_ins[1, 1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
                @assert isfinite(Σ_ins[1, 1, ki]) "fail Δ0 at $ki"
            end
        end
    end

    Σ.dynamic, Σ.instant = Σ_dyn ./ (-4 * π^2), Σ_ins ./ (-4 * π^2)
    return Σ
end

function G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
    @unpack dim = param
    # kernel = SelfEnergy.LegendreInteraction.DCKernel_old(param;
    # Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
    # kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
    #     Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
    # G0 = G0wrapped(Euv, rtol, kernel.kgrid, param)
    if dim == 2
        kernel = SelfEnergy.LegendreInteraction.DCKernel_2d(param;
            Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma)
        G0 = G0wrapped(Euv, rtol, kernel.kgrid, param)
        Σ = calcΣ_2d(G0, kernel)
    elseif dim == 3
        kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
            Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma)
        G0 = G0wrapped(Euv, rtol, kernel.kgrid, param)
        Σ = calcΣ_3d(G0, kernel)
    else
        error("No support for G0W0 in $dim dimension!")
    end

    return Σ
end

function G0W0(param; Euv=10 * param.EF, rtol=1e-14, Nk=12, maxK=6 * param.kF, minK=1e-8 * param.kF, order=4, int_type=:rpa)
    return G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
end

function zfactor(Σ::GreenFunc.Green2DLR)
    kgrid = Σ.spaceGrid
    kF = kgrid.panel[3]
    β = Σ.dlrGrid.β

    kF_label = searchsortedfirst(kgrid.grid, kF)
    Σ_freq = GreenFunc.toMatFreq(Σ, [0, 1])

    ΣI = imag(Σ_freq.dynamic[1, 1, kF_label, :])
    Z0 = 1 / (1 - (ΣI[2] - ΣI[1]) / 2 / π * β)
    # Z0 = 1 / (1 - imag(Σ_freq.dynamic[1, 1, kF_label, 1]) / π * β)

    return Z0
end

function massratio(param, Σ::GreenFunc.Green2DLR, δK=1e-6)
    # one can achieve ~1e-4 accuracy with δK = 1e-5
    @unpack kF, me = param

    δK *= kF
    kgrid = Σ.spaceGrid
    kF_label = searchsortedfirst(kgrid.grid, kF)
    z = zfactor(Σ)

    Σ_freq = GreenFunc.toMatFreq(Σ, [0, 1])
    k1, k2 = kF_label, kF_label + 1
    while abs(kgrid.grid[k2] - kgrid.grid[k1]) < δK
        k2 += 1
    end
    @assert kF < kgrid.grid[k1] < kgrid.grid[k2] "k1 and k2 are not on the same side! It breaks $kF > $(kgrid.grid[k1]) > $(kgrid.grid[k2])"
    sigma1 = real(Σ_freq.dynamic[1, 1, k1, 1] + Σ_freq.instant[1, 1, k1])
    sigma2 = real(Σ_freq.dynamic[1, 1, k2, 1] + Σ_freq.instant[1, 1, k2])
    ds_dk = (sigma1 - sigma2) / (kgrid.grid[k1] - kgrid.grid[k2])

    # println("m/kF ds_dk = $(me/kF*ds_dk)")
    return 1.0 / z / (1 + me / kF * ds_dk)
end

end
