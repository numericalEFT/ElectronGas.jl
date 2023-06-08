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
    @unpack me, β, μ = param

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol)
    green = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=ComplexF64)
    for ind in eachindex(green)
        green[ind] = 1 / (im * wn_mesh[ind[1]] - (green.mesh[2][ind[2]]^2 / 2 / me - μ))
    end

    return green
end

function Gwrapped(Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, param)
    @unpack me, kF, β, μ = param

    Σ_freq = Σ |> to_dlr |> to_imfreq
    green = similar(Σ_freq)

    w0i_label = locate(Σ_freq.mesh[1], 0)
    kf_label = locate(Σ_freq.mesh[2], kF)
    Σ_shift = real(Σ_freq[w0i_label, kf_label] + Σ_ins[1, kf_label])

    for ind in eachindex(green)
        green[ind] = 1 / (im * green.mesh[1][ind[1]] - (green.mesh[2][ind[2]]^2 / 2 / me - μ)
                          -
                          Σ_freq[ind] - Σ_ins[1, ind[2]] + Σ_shift)
    end

    return green
end

function calcΣ_2d(G::GreenFunc.MeshArray, W::LegendreInteraction.DCKernel)
    @unpack β = W.param

    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = G.mesh[1].representation
    bdlr = W.dlrGrid

    G_dlr = G |> to_dlr
    G_imt = G_dlr |> to_imtime

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3)

    # container of Σ
    Σ = GreenFunc.MeshArray(G_imt.mesh[1], kgrid; dtype=ComplexF64)

    # equal-time green (instant)
    G_ins = dlr_to_imtime(G_dlr, [β,]) * (-1)
    Σ_ins = GreenFunc.MeshArray(G_ins.mesh[1], kgrid; dtype=ComplexF64)

    for τi in eachindex(G_imt.mesh[1])
        for ki in eachindex(kgrid)
            Gq = CompositeGrids.Interp.interp1DGrid(G_imt[τi, :], G_imt.mesh[2], qgrids[ki].grid)
            integrand = kernel[ki, 1:qgrids[ki].size, τi] .* Gq .* qgrids[ki].grid
            Σ[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            @assert isfinite(Σ[τi, ki]) "fail Δ at $τi, $ki"
            if τi == 1
                Gq = CompositeGrids.Interp.interp1DGrid(G_ins[1, :], G_ins.mesh[2], qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Gq .* qgrids[ki].grid
                Σ_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
                @assert isfinite(Σ_ins[1, ki]) "fail Δ0 at $ki"
            end
        end
    end

    return Σ / (-4 * π^2), Σ_ins / (-4 * π^2)
end

function calcΣ_3d(G::GreenFunc.MeshArray, W::LegendreInteraction.DCKernel)
    @unpack β = W.param

    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = G.mesh[1].representation
    bdlr = W.dlrGrid

    G_dlr = G |> to_dlr
    G_imt = G_dlr |> to_imtime

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3)

    # container of Σ
    Σ = GreenFunc.MeshArray(G_imt.mesh[1], kgrid; dtype=ComplexF64)

    # equal-time green (instant)
    G_ins = dlr_to_imtime(G_dlr, [β,]) * (-1)
    Σ_ins = GreenFunc.MeshArray(G_ins.mesh[1], kgrid; dtype=ComplexF64)

    for τi in eachindex(G_imt.mesh[1])
        for ki in eachindex(kgrid)
            k = kgrid[ki]
            Gq = CompositeGrids.Interp.interp1DGrid(G_imt[τi, :], G_imt.mesh[2], qgrids[ki].grid)
            integrand = kernel[ki, 1:qgrids[ki].size, τi] .* Gq ./ k .* qgrids[ki].grid
            Σ[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            @assert isfinite(Σ[τi, ki]) "fail Δ at $τi, $ki"
            if τi == 1
                Gq = CompositeGrids.Interp.interp1DGrid(G_ins[1, :], G_ins.mesh[2], qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Gq ./ k .* qgrids[ki].grid
                Σ_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
                @assert isfinite(Σ_ins[1, ki]) "fail Δ0 at $ki"
            end
        end
    end

    return Σ / (-4 * π^2), Σ_ins / (-4 * π^2)
end

function G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type, kgrid::Union{AbstractGrid,AbstractVector,Nothing}=nothing; kwargs...)
    @unpack dim = param
    # kernel = SelfEnergy.LegendreInteraction.DCKernel_old(param;
    # Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
    # kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
    #     Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
    # G0 = G0wrapped(Euv, rtol, kernel.kgrid, param)

    kGgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)
    if dim == 2
        if isnothing(kgrid)
            kernel = SelfEnergy.LegendreInteraction.DCKernel_2d(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma, kwargs...)
        else
            if (kgrid isa AbstractVector)
                kgrid = SimpleG.Arbitrary{eltype(kgrid)}(kgrid)
            end
            kernel = SelfEnergy.LegendreInteraction.DCKernel_2d(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma, kgrid=kgrid, kwargs...)
        end
        G0 = G0wrapped(Euv, rtol, kGgrid, param)
        Σ, Σ_ins = calcΣ_2d(G0, kernel)
    elseif dim == 3
        if isnothing(kgrid)
            kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma, kwargs...)
        else
            if (kgrid isa AbstractVector)
                kgrid = SimpleG.Arbitrary{eltype(kgrid)}(kgrid)
            end
            kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma, kgrid=kgrid, kwargs...)
        end
        G0 = G0wrapped(Euv, rtol, kGgrid, param)
        Σ, Σ_ins = calcΣ_3d(G0, kernel)
    else
        error("No support for G0W0 in $dim dimension!")
    end

    return Σ, Σ_ins
end

function G0W0(param, kgrid::Union{AbstractGrid,AbstractVector,Nothing}=nothing; Euv=100 * param.EF, rtol=1e-14, Nk=12, maxK=6 * param.kF, minK=1e-8 * param.kF, order=8, int_type=:rpa,
    kwargs...)
    return G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type, kgrid; kwargs...)
end

function G0Γ0(param, kgrid::Union{AbstractGrid,AbstractVector,Nothing}=nothing; Euv=100 * param.EF, rtol=1e-14, Nk=12, maxK=6 * param.kF, minK=1e-8 * param.kF, order=8, int_type=:none,
    kwargs...)
    return G0Γ0(param, Euv, rtol, Nk, maxK, minK, order, int_type, kgrid; kwargs...)
end

function G0Γ0(param, Euv, rtol, Nk, maxK, minK, order, int_type, kgrid::Union{AbstractGrid,AbstractVector,Nothing}=nothing; kwargs...)
    @unpack dim = param
    # kernel = SelfEnergy.LegendreInteraction.DCKernel_old(param;
    # Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
    # kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
    #     Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
    # G0 = G0wrapped(Euv, rtol, kernel.kgrid, param)

    kGgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)
    if dim == 2
        error("No support for G0Γ0 in 2d dimension!")
    elseif dim == 3
        if isnothing(kgrid)
            kernel = SelfEnergy.LegendreInteraction.DCKernel_Ladder(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma, kwargs...)
        else
            if (kgrid isa AbstractVector)
                kgrid = SimpleG.Arbitrary{eltype(kgrid)}(kgrid)
            end
            kernel = SelfEnergy.LegendreInteraction.DCKernel_Ladder(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma, kgrid=kgrid, kwargs...)
        end
        G0 = G0wrapped(Euv, rtol, kGgrid, param)
        Σ, Σ_ins = calcΣ_3d(G0, kernel)
    else
        error("No support for G0W0 in $dim dimension!")
    end

    return Σ
end

"""
    function zfactor(param, Σ::GreenFunc.MeshArray; kamp=param.kF, ngrid=[0, 1])
    
calculate the z-factor of the self-energy at the momentum kamp
```math
    z_k=\\frac{1}{1-\\frac{\\partial Im\\Sigma(k, 0^+)}{\\partial \\omega}}
```
"""
function zfactor(param, Σ::GreenFunc.MeshArray; kamp=param.kF, ngrid=[0, 1])
    @unpack kF, β = param

    k_label = locate(Σ.mesh[2], kamp)
    kamp = Σ.mesh[2][k_label]

    Σ_freq = dlr_to_imfreq(to_dlr(Σ), ngrid[1:2])
    ΣI = imag(Σ_freq[:, k_label])
    ds_dw = (ΣI[2] - ΣI[1]) / 2 / π * β
    Z0 = 1 / (1 - ds_dw)

    return Z0, kamp
end

"""
    function massratio(param, Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, δK=5e-6; kamp=param.kF)
    
calculate the effective mass of the self-energy at the momentum kamp
```math
    \\frac{m^*_k}{m}=\\frac{1}{z_k} \\cdot \\left(1+\\frac{m}{k}\\frac{\\partial Re\\Sigma(k, 0)}{\\partial k}\\right)^{-1}
```
"""
function massratio(param, Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, δK=5e-6; kamp=param.kF)
    # one can achieve ~1e-5 accuracy with δK = 5e-6
    @unpack kF, me = param

    δK *= kF
    k_label = locate(Σ.mesh[2], kamp)
    kamp = Σ.mesh[2][k_label]
    z = zfactor(param, Σ; kamp=kamp)[1]

    Σ_freq = dlr_to_imfreq(to_dlr(Σ), [0, 1])
    k1, k2 = k_label, k_label + 1
    while abs(Σ.mesh[2][k2] - Σ.mesh[2][k1]) < δK
        k2 += 1
    end
    # @assert kF < kgrid.grid[k1] < kgrid.grid[k2] "k1 and k2 are not on the same side! It breaks $kF > $(kgrid.grid[k1]) > $(kgrid.grid[k2])"
    sigma1 = real(Σ_freq[1, k1] + Σ_ins[1, k1])
    sigma2 = real(Σ_freq[1, k2] + Σ_ins[1, k2])
    ds_dk = (sigma1 - sigma2) / (Σ.mesh[2][k1] - Σ.mesh[2][k2])

    return 1.0 / z / (1 + me / kamp * ds_dk), kamp
end

"""
    function bandmassratio(param, Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray; kamp=param.kF)
    
calculate the effective band mass of the self-energy at the momentum kamp
```math
    \\frac{m^*_k}{m}=\\frac{1}{z_k}\\cdot \\left(1+\\frac{Re\\Sigma(k, 0) - Re\\Sigma(0, 0)}{k^2/2m}\\right)^{-1}
```
"""
function bandmassratio(param, Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray; kamp=param.kF)
    # one can achieve ~1e-5 accuracy with δK = 5e-6
    @unpack me = param
    z = zfactor(param, Σ; kamp=kamp)[1]

    k_label = locate(Σ.mesh[2], kamp)
    kamp = Σ.mesh[2][k_label]

    Σ_freq = dlr_to_imfreq(to_dlr(Σ), [0, 1])
    sigma1 = real(Σ_freq[1, k_label] + Σ_ins[1, k_label])
    sigma2 = real(Σ_freq[1, 1] + Σ_ins[1, 1])
    ds_dk = (sigma1 - sigma2) / (kamp^2 / 2 / me)

    return 1.0 / z / (1 + ds_dk), kamp
end

function chemicalpotential(param, Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray)
    # one can achieve ~1e-5 accuracy with δK = 5e-6
    @unpack kF, me = param
    k_label = locate(Σ.mesh[2], kF)

    Σ_freq = dlr_to_imfreq(to_dlr(Σ), [-1, 0])
    return real(Σ_freq[1, k_label] + Σ_ins[1, k_label])
end

end
