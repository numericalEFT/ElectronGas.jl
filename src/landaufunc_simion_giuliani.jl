
abstract type AbstractFitCoeffs end

# Perdew-Wang fit coefficients for α_c(rs), α^{RPA}_c(rs), E_c(rs), and E^{RPA}_c(rs).
# See Table I of J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992) (doi.org/10.1103/PhysRevB.45.13244).
struct PerdewWangFitCoeffs <: AbstractFitCoeffs
    p::Float64
    Ac::Float64
    α1::Float64
    β1::Float64
    β2::Float64
    β3::Float64
    β4::Float64
end

# Perdew-Wang interpolating function for α_c(rs), α^{RPA}_c(rs), E_c(rs), and E^{RPA}_c(rs)
@inline function perdew_wang_interpolant(fitcoeffs::PerdewWangFitCoeffs)
    @unpack p, Ac, α1, β1, β2, β3, β4 = fitcoeffs
    interpolant(rs) = -4Ac * (1 + α1 * rs) * log1p(1 / (2Ac * (β1 * sqrt(rs) + β2 * rs + β3 * rs^(3 / 2.0) + β4 * rs^(p + 1.0))))
    # interpolant(rs) = -4Ac * (1 + α1 * rs) * log(1 + 1 / (2Ac * (β1 * sqrt(rs) + β2 * rs + β3 * rs^(3 / 2.0) + β4 * rs^(p + 1.0))))
    return interpolant
end

# Perdew-Wang fit to E^{RPA}_c(rs) for ζ = 0 (paramagnetic phase)
@inline const E_corr_rpa_z0 = perdew_wang_interpolant(
    PerdewWangFitCoeffs(0.75, 0.031091, 0.082477, 5.1486, 1.6483, 0.23647, 0.20614)
)

# Perdew-Wang fit to E^{RPA}_c(rs) for ζ = 1 (ferromagnetic phase)
@inline const E_corr_rpa_z1 = perdew_wang_interpolant(
    PerdewWangFitCoeffs(0.75, 0.015545, 0.035374, 6.4869, 1.3083, 0.1518, 0.082349)
)

# Perdew-Wang fit to -α^{RPA}_c(rs)
@inline const negative_spin_stiffness_rpa = perdew_wang_interpolant(
    PerdewWangFitCoeffs(1.0, 0.016887, 0.028829, 10.357, 3.6231, 0.4799, 0.12279)
)

# Perdew-Wang fit to E_c(rs) for ζ = 0 (paramagnetic phase)
@inline const E_corr_z0 = perdew_wang_interpolant(
    PerdewWangFitCoeffs(1.0, 0.031091, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294)
)

# Perdew-Wang fit to E_c(rs) for ζ = 1 (ferromagnetic phase)
@inline const E_corr_z1 = perdew_wang_interpolant(
    PerdewWangFitCoeffs(1.0, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
)

# Perdew-Wang fit to -α_c(rs)
@inline const negative_spin_stiffness = perdew_wang_interpolant(
    PerdewWangFitCoeffs(1.0, 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)
)

@inline function spin_susceptibility_enhancement(rs)
    alpha_ueg = (4 / 9π)^(1 / 3)
    # return 1 - (alpha_ueg * rs / π) + 3(alpha_ueg * rs)^2 * spin_stiffness(rs)
    # return 1 - (alpha_ueg * rs / π) - 3(alpha_ueg * rs)^2 * negative_spin_stiffness(rs)
    return 1 - (alpha_ueg * rs / π) - 3(alpha_ueg * rs)^2 * negative_spin_stiffness(rs) / 2
end

# F_s from Moroni (doi:10.1103/PhysRevB.57.14569)
@inline function F_s_Moroni(q, param::Parameter.Para)
    @unpack kF, rs, e0 = param
    n0 = (kF)^3 / 3 / π^2
    if e0 ≈ 0.0
        return 0.0
    end

    # step for numerical derivatives 
    step = rs < 1.0 ? 1e-3 : 1e-7  # NOTE: Numerically unstable near rs = 0! TODO: Fix the instability
    # step = 1e-7
    x = sqrt(rs)

    # Calculate parameter A
    rs1 = (n0 / (n0 + step))^(1.0 / 3.0) * rs
    rs_1 = (n0 / (n0 - step))^(1.0 / 3.0) * rs
    deriv_1 = (E_corr_z0(rs1) - E_corr_z0(rs)) / step
    deriv_2 = (E_corr_z0(rs1) + E_corr_z0(rs_1) - 2 * E_corr_z0(rs)) / step^2
    A = 0.25 - kF^2 / 4 / π / e0^2 * (2 * deriv_1 + n0 * deriv_2)
    #println("A=$(A)")

    # Calculate parameter B
    a1 = 2.15
    a2 = 0.435
    b1 = 1.57
    b2 = 0.409
    B = (1 + a1 * x + a2 * x^3) / (3 + b1 * x + b2 * x^3)
    #println("B=$(B)")

    # Calculate parameter C
    step = 1e-7  # NOTE: no instability in deriv_1 near rs = 0!
    deriv_1 = (E_corr_z0(rs + step) - E_corr_z0(rs)) / step
    #deriv_1 = E_corr_p(rs)
    C = -π / 2 / kF / e0^2 * (E_corr_z0(rs) + rs * deriv_1)
    #println("C=$(C)")

    D = B / (A - C)
    α = 1.5 / rs^0.25 * A / B / D
    β_0 = 1.2 / B / D
    Q = q / kF
    G_s = C * Q^2 + B * Q^2 / (D + Q^2) + α * Q^4 * exp(-β_0 * Q^2)
    F_s = 4 * π * e0^2 * G_s / q^2
    return F_s
end

# F_a from Simion & Giuliani (doi.org/10.1103/PhysRevB.77.035131)
@inline function F_a_Simion_Giuliani(q, param::Parameter.Para)
    @unpack kF, qTF, rs, e0 = param
    if e0 ≈ 0.0
        return 0.0
    end

    # step for numerical derivatives 
    step = 1e-7  # NOTE: no instability in deriv_1 near rs = 0!
    x = sqrt(rs)

    # Calculate parameter B1
    a1 = 2.15
    a2 = 0.435
    b1 = 1.57
    b2 = 0.409
    B = (1 + a1 * x + a2 * x^3) / (3 + b1 * x + b2 * x^3)
    g0 = 32.0 / (8.0 + 3.0 * rs)^2
    B1 = B - 1 + 2g0  # Eq. (34)
    #println("B1=$(B1)")

    # Calculate parameter C
    deriv_1 = (E_corr_z0(rs + step) - E_corr_z0(rs)) / step
    #deriv_1 = E_corr_p(rs)
    C = -π / 2 / kF / e0^2 * (E_corr_z0(rs) + rs * deriv_1)
    #println("C=$(C)")

    # Calculate parameter D1
    chi0_over_chi = spin_susceptibility_enhancement(rs)
    D1 = B1 / ((kF^2 / qTF^2) * (1.0 - chi0_over_chi) - C)  # Corrected Eq. (36): κ² = 4πe²N(0) = qTF²
    #println("D1=$(D1)")

    Q = q / kF
    G_a = C * Q^2 + B1 * Q^2 / (D1 + Q^2)  # Eq. (33)
    F_a = 4 * π * e0^2 * G_a / q^2  # F(+-) = G(+-) * V
    return F_a
end

"""
    function landauParameterSimionGiuliani(q, n, param; kwargs...)

Ansatz for spin symmetric and antisymmetric Landau parameters F_s and F_a
following Simion & Giuliani (doi.org/10.1103/PhysRevB.77.035131).

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function landauParameterSimionGiuliani(q, n, param; kwargs...)
    F_s = F_s_Moroni(q, param)
    F_a = F_a_Simion_Giuliani(q, param)
    return F_s, F_a
end
