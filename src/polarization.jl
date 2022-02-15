module Polarization

using ..GreenFunc, ..CompositeGrids
using ..Parameter
using ..Parameters
using ..Convention

export Polarization0_ZeroTemp, Polarization0_FiniteTemp

# Analytical calculated integrand of Π0 in 2D.
@inline function _ΠT2d_integrand(k, q, ω, param)
    @unpack me, β, μ = param
    density = me / 2π
    nk = 1.0 / (exp(β * (k^2 / 2 / me - μ)) + 1)

    return nothing
end

# Analytical calculated integrand of Π0 in 3D.
@inline function _ΠT3d_integrand(k, q, ω, param)
    @unpack me, β, μ = param
    # ω only appears as ω^2 so no need to check sign of ω
    nk = 1.0 / (exp(β * (k^2 / 2 / me - μ)) + 1)

    # if q is too small, use safe form
    if q < 1e-16 && ω ≈ 0
        if abs(q - 2 * k)^2 < 1e-16
            return 0.0
        else
            return -k * me / (4 * π^2) * nk * ((8 * k) / ((q - 2 * k)^2))
        end
    elseif q < 1e-16 && !(ω ≈ 0)
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

function finitetemp_kgrid(q::Float64, kF::Float64, maxk = 20, scaleN = 20, minterval = 1e-6, gaussN = 10)
    mink = (q < 1e-16 / minterval) ? minterval * kF : minterval * min(q, kF)
    kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxk * kF], [0.5 * q, kF], scaleN, mink, gaussN)
    return kgrid
end

"""
    function Polarization0_FiniteTemp(q, n, param, maxk=20, scaleN=20, minterval=1e-6, gaussN=10)

Finite temperature one-spin Π0 function for matsubara frequency and momentum. Analytically sum over transfer frequency and angular
dependence of momentum, and numerically calculate integration of magnitude of momentum.
Slower(~200μs) than Polarization0_ZeroTemp.
Assume G_0^{-1} = iω_n - (k^2/(2m) - mu)

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
 - maxk: optional, upper limit of integral -> maxk*kF
 - scaleN: optional, N of Log grid in LogDensedGrid, check CompositeGrids for more detail
 - minterval: optional, actual minterval of grid is this value times min(q,kF)
 - gaussN: optional, N of GaussLegendre grid in LogDensedGrid.
"""
function Polarization0_FiniteTemp(q::Float64, n::Int, param; maxk = 20, scaleN = 20, minterval = 1e-6, gaussN = 10)
    @unpack dim, kF, β = param
    if dim ∉ [2, 3]
        error("No support for finite-temperature polarization in $dim dimension!")
    end
    # check sign of q, use -q if negative
    if q < 0
        q = -q
    end

    # mink = (q < 1e-16 / minterval) ? minterval * kF : minterval * min(q, kF)
    # kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxk * kF], [0.5 * q, kF], scaleN, mink, gaussN)
    kgrid = finitetemp_kgrid(q, kF, maxk, scaleN, minterval, gaussN)
    integrand = zeros(Float64, kgrid.size)
    if dim == 2
        for (ki, k) in enumerate(kgrid.grid)
            integrand[ki] = _ΠT2d_integrand(k, q, 2π * n / β, param)
            @assert !isnan(integrand[ki]) "nan at k=$k, q=$q"
        end
    elseif dim == 3
        for (ki, k) in enumerate(kgrid.grid)
            integrand[ki] = _ΠT3d_integrand(k, q, 2π * n / β, param)
            @assert !isnan(integrand[ki]) "nan at k=$k, q=$q"
        end
    end

    return Interp.integrate1D(integrand, kgrid)
end

function Polarization0_FiniteTemp(q::Float64, n::AbstractVector, param; maxk = 20, scaleN = 20, minterval = 1e-6, gaussN = 10)
    @unpack dim, kF, β = param
    if dim ∉ [2, 3]
        error("No support for finite-temperature polarization in $dim dimension!")
    end
    # check sign of q, use -q if negative
    if q < 0
        q = -q
    end

    # mink = (q < 1e-16 / minterval) ? minterval * kF : minterval * min(q, kF)
    # kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxk * kF], [0.5 * q, kF], scaleN, mink, gaussN)
    kgrid = finitetemp_kgrid(q, kF, maxk, scaleN, minterval, gaussN)
    integrand = zeros(Float64, (kgrid.size, length(n)))
    if dim == 2
        for (ki, k) in enumerate(kgrid.grid)
            for (mi, m) in enumerate(n)
                integrand[ki, mi] = _ΠT2d_integrand(k, q, 2π * m / β, param)
                @assert !isnan(integrand[ki, mi]) "nan at k=$k, q=$q"
            end
        end
    elseif dim == 3
        for (ki, k) in enumerate(kgrid.grid)
            for (mi, m) in enumerate(n)
                integrand[ki, mi] = _ΠT3d_integrand(k, q, 2π * m / β, param)
                @assert !isnan(integrand[ki, mi]) "nan at k=$k, q=$q"
            end
        end
    end

    return Interp.integrate1D(integrand, kgrid; axis = 1)
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
    ω_n = 2 * π * n / β
    v = me * ω_n / q / kF * im + q / 2 / kF
    # v2 = me * ω_n / q / kF * im - q / 2 / kF

    if q < EPS && n != 0
        Π = 0.0
    else
        Π = 1.0 - real(√(v^2 - 1)) * 2kF / q
        # Π = 1.0 - real(√(v^2 - 1) + √(v2^2 - 1)) * kF / q
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
    if n < 0
        n = -n
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

function Polarization0_ZeroTemp(q::Float64, n::AbstractVector, param)
    @unpack dim = param
    result = zeros(Float64, length(n))
    if dim == 2
        for (mi, m) in enumerate(n)
            result[mi] = Polarization0_2dZeroTemp(q, m, param)
        end
    elseif dim == 3
        for (mi, m) in enumerate(n)
            result[mi] = Polarization0_3dZeroTemp(q, m, param)
        end
    else
        error("No support for zero-temperature polarization in $dim dimension!")
    end
    return result
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
 - polatype: type of pi function, support :zerotemp or :finitetemp
"""
function Polarization0wrapped(Euv, rtol, sgrid::SGT, param; pifunc = Polarization0_ZeroTemp) where {SGT}
    @unpack β = param

    green = GreenFunc.Green2DLR{Float64}(:polarization, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry = :ph, rtol = rtol)
    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(sgrid)
        # for (ni, n) in enumerate(green.dlrGrid.n)
        #     green_dyn[1, 1, ki, ni] = pifunc(k, n, param)
        # end
        green_dyn[1, 1, ki, :] = pifunc(k, green.dlrGrid.n, param)
    end
    green.dynamic = green_dyn
    return green
end

end
