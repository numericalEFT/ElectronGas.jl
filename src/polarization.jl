module Polarization

using ..GreenFunc, ..CompositeGrids
using ..Parameter
using ..Parameters
using ..Convention

export Polarization0_ZeroTemp, Polarization0_FiniteTemp

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd() * "/" * ARGS[1])

# using Parameters, GreenFunc, CompositeGrids
# include(srcdir*"/parameter.jl")
# using .Parameter
# include(srcdir*"/convention.jl")
# using .Convention

# Analytical calculated integrand of Π0 in 2D.
# @inline function _ΠT2d_integrand(k, q, ω, param)
#     @unpack me, β, EF = param
#     density = me / 2π

#     v1 = me * ω - q^2 / 2
#     v2 = me * ω + q^2 / 2
#     nk = 1.0 / (exp(β * (k^2 / 2 / me - EF)) + 1)

#     theta1, theta2 = abs(v1) - k * q, abs(v2) - k * q
#     if theta1 <= 0 && theta2 <= 0
#         integrand = 0.0
#     elseif theta1 <= 0 && theta2 > 0
#         integrand = k * nk * (-sign(v2) / √(v2^2 - k^2 * q^2))
#     elseif theta1 > 0 && theta2 <= 0
#         integrand = k * nk * sign(v1) / √(v1^2 - k^2 * q^2)
#     else
#         integrand = k * nk * (sign(v1) / √(v1^2 - k^2 * q^2) - sign(v2) / √(v2^2 - k^2 * q^2))
#     end

#     return density * integrand
# end

# Analytical calculated integrand of Π0 in 2D.
@inline function _ΠT2d_integrand(k, q, ω, param)
    @unpack me, β, EF = param
    density = me / 2π
    nk = 1.0 / (exp(β * (k^2 / 2 / me - EF)) + 1)

    # if q is too small, use safe form
    if q < EPS && ω == 0
        if abs(q - 2 * k)^2 < EPS
            return 0.0
        else
            return -density / 2 * nk * nk * ((8 * k) / ((q - 2 * k)^2))
        end
    elseif q < EPS && ω != 0
        return -density / 2 * nk * ((8 * k * q^2) / (4 * me^2 * ω^2 + (q^2 - 2 * k * q)^2))
    else
        return -density / 2q * nk * log1p((8 * k * q^3) / (4 * me^2 * ω^2 + (q^2 - 2 * k * q)^2))
    end

end

# Analytical calculated integrand of Π0 in 3D.
@inline function _ΠT3d_integrand(k, q, ω, param)
    @unpack me, β, EF = param
    # ω only appears as ω^2 so no need to check sign of ω
    nk = 1.0 / (exp(β * (k^2 / 2 / me - EF)) + 1)

    # if q is too small, use safe form
    if q < 1e-16 && ω == 0
        if abs(q - 2 * k)^2 < 1e-16
            return 0.0
        else
            return -k * me / (4 * π^2) * nk * ((8 * k) / ((q - 2 * k)^2))
        end
    elseif q < 1e-16 && ω != 0
        return -k * me / (4 * π^2) * nk * ((8 * k * q^2) / (4 * me^2 * ω^2 + (q^2 - 2 * k * q)^2))
    else
        return -k * me / (4 * π^2 * q) * nk * log1p((8 * k * q^3) / (4 * me^2 * ω^2 + (q^2 - 2 * k * q)^2))
    end
    # if ω==0 && abs(q*(k*2-q)*β)<1e-16
    #     return 0.0
    # elseif ω!=0 && k^2*q^2/4/me^2<ω^2*1e-16
    #     2*q^2*k^2/(π^2*ω^2)/(exp(β*(k^2/2/me-μ))+1)/(8*me)
    # else
    #     return k*me/(2*π^2*q)/(exp(β*(k^2/2/me-μ))+1)*log((4*me^2*ω^2+(q^2+2*k*q)^2)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    #     return k*me/(2*π^2*q)/(exp(β*(k^2/2/me-μ))+1)*log1p((8*k*q^3)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    # end
end

"""
    function Polarization0_FiniteTemp(q, n, param, maxk=20, scaleN=20, minterval=1e-6, gaussN=10)

Finite temperature one-spin Π0 function for matsubara frequency and momentum. Analytically sum over transfer frequency and angular
dependence of momentum, and numerically calculate integration of magnitude of momentum.
Slower(~200μs) than Polarization0_ZeroTemp.
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F)

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
 - maxk: optional, upper limit of integral -> maxk*kF
 - scaleN: optional, N of Log grid in LogDensedGrid, check CompositeGrids for more detail
 - minterval: optional, actual minterval of grid is this value times min(q,kF)
 - gaussN: optional, N of GaussLegendre grid in LogDensedGrid.
"""
function Polarization0_FiniteTemp(q::Float64, n::Int, param, maxk = 20, scaleN = 20, minterval = 1e-6, gaussN = 10)
    @unpack dim, kF, β = param
    if dim ∉ [2, 3]
        error("No support for finite-temperature polarization in $dim dimension!")
    end
    # check sign of q, use -q if negative
    if q < 0
        q = -q
    end

    if dim == 2
        mink = (q < 1e-16 / minterval) ? minterval * kF : minterval * min(q, kF)
        kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxk * kF], [0.5 * q, kF], scaleN, mink, gaussN)
        integrand = zeros(Float64, kgrid.size)
        for (ki, k) in enumerate(kgrid.grid)
            integrand[ki] = _ΠT2d_integrand(k, q, 2π * n / β, param)
            @assert !isnan(integrand[ki]) "nan at k=$k, q=$q"
        end
    elseif dim == 3
        mink = (q < 1e-16 / minterval) ? minterval * kF : minterval * min(q, kF)
        kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxk * kF], [0.5 * q, kF], scaleN, mink, gaussN)
        integrand = zeros(Float64, kgrid.size)
        for (ki, k) in enumerate(kgrid.grid)
            integrand[ki] = _ΠT3d_integrand(k, q, 2π * n / β, param)
            @assert !isnan(integrand[ki]) "nan at k=$k, q=$q"
        end
    end

    return Interp.integrate1D(integrand, kgrid)
end


@inline function Polarization0_2dZeroTemp(q, n, param)
    @unpack me, kF, β = param
    density = me / 2π
    # check sign of q, use -q if negative
    if q < 0
        q = -q
    end
    # if q is too small, set to 1000eps
    if q < eps(0.0) * 1e6
        q = eps(0.0) * 1e6
    end

    Π = 0.0
    x = q / 2 / kF
    ω_n = 2 * π * n / β
    y = me * ω_n / q / kF

    if abs(y - x) <= 1 && abs(y + x) <= 1
        Π = 1.0
    elseif abs(y - x) <= 1 && abs(y + x) > 1
        Π = 1.0 - sign(y + x) * √((y + x)^2 - 1) / 2x
    elseif abs(y - x) > 1 && abs(y + x) <= 1
        Π = 1.0 + sign(y - x) * √((y - x)^2 - 1) / 2x
    else
        z = q / (me * ω_n)
        a = 1.0 / (me * ω_n)
        if z < 3.8e-4
            a = 1.0 / (me * ω_n)
            Π = -z^2 / 2 - 3 * z^4 / 8 - (5 + 8 * a^2) * z^6 / 16
        else
            Π = 1.0 + (sign(y - x) * √((y - x)^2 - 1) - sign(y + x) * √((y + x)^2 - 1)) / 2x
        end
    end

    return -density * Π
end

@inline function Polarization0_3dZeroTemp(q, n, param)
    @unpack me, kF, β = param
    density = me * kF / (2π^2)
    # check sign of q, use -q if negative
    if q < 0
        q = -q
    end
    # if q is too small, set to 1000eps
    if q < eps(0.0) * 1e6
        q = eps(0.0) * 1e6
    end

    Π = 0.0
    x = q / 2 / kF
    ω_n = 2 * π * n / β
    y = me * ω_n / q / kF

    if n == 0
        if abs(q - 2 * kF) > EPS
            Π = density * (1 + (1 - x^2) * log1p(4 * x / ((1 - x)^2)) / 4 / x)
        else
            Π = density
        end
    else
        if abs(q - 2 * kF) > EPS
            if y^2 < 1e-4 / EPS
                theta = atan(2 * y / (y^2 + x^2 - 1))
                if theta < 0
                    theta = theta + π
                end
                @assert theta >= 0 && theta <= π
                Π = density * (1 + (1 - x^2 + y^2) * log1p(4 * x / ((1 - x)^2 + y^2)) / 4 / x - y * theta)
            else
                Π = density * (2.0 / 3.0 / y^2 - 2.0 / 5.0 / y^4) #+ (6.0 - 14.0*(ω_n/4.0)^2)/21.0/y^6)
            end
        else
            theta = atan(2 / y)
            if theta < 0
                theta = theta + π
            end
            @assert theta >= 0 && theta <= π
            Π = density * (1 + y^2 * log1p(4 / (y^2)) / 4 - y * theta)
        end
    end
    # initially derived for spin=1/2
    return -Π / 2
end

"""
    function Polarization0_ZeroTemp(q, n, param)

Zero temperature one-spin Π0 function for matsubara frequency and momentum. For low temperature the finite temperature
polarization could be approximated with this function to run faster(~200ns).
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F).


#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
  - param: other system parameters
"""
function Polarization0_ZeroTemp(q::Float64, n::Int, param)
    @unpack dim = param

    if dim == 2
        return Polarization0_2dZeroTemp(q, n, param)
    elseif dim == 3
        return Polarization0_3dZeroTemp(q, n, param)
    else
        error("No support for zero-temperature polarization in $dim dimension!")
    end
end


"""
    function Polarization0wrapped(Euv, rtol, sgrid::SGT, param, polatype=:zerotemp) where{TGT, SGT}

Π0 function for matsubara frequency and momentum. Use Polarization0_ZeroTemp by default,
Polarization0_FiniteTemp when pifunc is specified.
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F).
Return full polarization0 function stored in GreenFunc.GreenBasic.Green2DLR.

#Arguments:
 - Euv: Euv of DLRGrid
 - rtol: rtol of DLRGrid
 - sgrid: momentum grid
 - param: other system parameters
 - pifunc: single point Π0 function used. Require form with pifunc(k, n, param).
"""
function Polarization0wrapped(Euv, rtol, sgrid::SGT, param, pifunc = Polarization0_ZeroTemp) where {SGT}
    @unpack β = param

    green = GreenFunc.Green2DLR{Float64}(:polarization, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry = :ph, rtol = rtol)
    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1, 1, ki, ni] = pifunc(k, n, param)
        end
    end
    green.dynamic = green_dyn
    return green
end

end
