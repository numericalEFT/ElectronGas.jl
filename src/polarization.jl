module Polarization

using Parameters, GreenFunc, CompositeGrids

export Polarization0_ZeroTemp, Polarization0_FiniteTemp

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

include(srcdir*"/parameter.jl")
using .Parameter

include(srcdir*"/convention.jl")
using .Convention

# Analytical calculated integrand of Π0.
@inline function _ΠT_integrand(k, q, ω, param)
    @unpack me, beta, EF, kF = param
    # ω only appears as ω^2 so no need to check sign of ω

    # if q is too small, use safe form
    if q < 1e-16 && ω==0
        if abs(q-2*k)^2<1e-16
            return 0.0
        else
            return k*me/(2*π^2)/(exp(beta*(k^2/2/me-EF))+1)*((8*k)/((q-2*k)^2))
        end
    elseif q < 1e-16 && ω!=0
        return k*me/(2*π^2)/(exp(beta*(k^2/2/me-EF))+1)*((8*k*q^2)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    else
        return k*me/(2*π^2*q)/(exp(beta*(k^2/2/me-EF))+1)*log1p((8*k*q^3)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    end
    # if ω==0 && abs(q*(k*2-q)*beta)<1e-16
    #     return 0.0
    # elseif ω!=0 && k^2*q^2/4/me^2<ω^2*1e-16
    #     2*q^2*k^2/(π^2*ω^2)/(exp(beta*(k^2/2/me-μ))+1)/(8*me)
    # else
    #     return k*me/(2*π^2*q)/(exp(beta*(k^2/2/me-μ))+1)*log((4*me^2*ω^2+(q^2+2*k*q)^2)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    #     return k*me/(2*π^2*q)/(exp(beta*(k^2/2/me-μ))+1)*log1p((8*k*q^3)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    # end
end

"""
    function Polarization0_FiniteTemp(q, n, param, maxk=20, scaleN=20, minterval=1e-6, gaussN=10)

Finite temperature Π0 function for matsubara frequency and momentum. Analytically sum over transfer frequency and angular
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
function Polarization0_FiniteTemp(q, n, param, maxk=20, scaleN=20, minterval=1e-6, gaussN=10)
    @unpack me, kF, beta = param
    # check sign of q, use -q if negative
    if q<0
        q = -q
    end
    mink = (q<1e-16/minterval) ? minterval*kF : minterval*min(q,kF)
    kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxk*kF], [0.5*q, kF], scaleN, mink, gaussN)
    integrand = zeros(Float64, kgrid.size)
    for (ki, k) in enumerate(kgrid.grid)
        integrand[ki] = _ΠT_integrand(k, q, 2π*n/beta, param)
        @assert !isnan(integrand[ki]) "nan at k=$k, q=$q"
    end

    return Interp.integrate1D(integrand, kgrid)
end

"""
    function Polarization0_ZeroTemp(q, n, param)

Zero temperature Π0 function for matsubara frequency and momentum. For low temperature the finite temperature
polarization could be approximated with this function to run faster(~200ns).
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F).


#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function Polarization0_ZeroTemp(q, n, param)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0 = param
    # check sign of q, use -q if negative
    if q<0
        q = -q
    end

    # if q is too small, set to 1000eps
    if q < eps(0.0)*1e6
        q = eps(0.0)*1e6
    end

    Π = 0.0
    x = q/2/kF
    ω_n = 2*π*n/beta
    y = me*ω_n/q/kF

    if n == 0
        if abs(q - 2*kF) > EPS
            Π = me*kF/2/π^2*(1 + (1 -x^2)*log1p(4*x/((1-x)^2))/4/x)
        else
            Π = me*kF/2/π^2
        end
    else
        if abs(q - 2*kF) > EPS
            if y^2 < 1e-4/EPS                    
                theta = atan( 2*y/(y^2+x^2-1) )
                if theta < 0
                    theta = theta + π
                end
                @assert theta >= 0 && theta<= π
                Π = me*kF/2/π^2 * (1 + (1 -x^2 + y^2)*log1p(4*x/((1-x)^2+y^2))/4/x - y*theta)
            else
                Π = me*kF/2/π^2 * (2.0/3.0/y^2  - 2.0/5.0/y^4) #+ (6.0 - 14.0*(ω_n/4.0)^2)/21.0/y^6)
            end
        else
            theta = atan( 2/y )
            if theta < 0
                theta = theta + π
            end
            @assert theta >= 0 && theta<= π
            Π = me*kF/2/π^2*(1 + y^2*log1p(4/(y^2))/4 - y*theta)
        end
    end

    return Π
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
function Polarization0wrapped(Euv, rtol, sgrid::SGT, param, pifunc=Polarization0_ZeroTemp) where{SGT}
    @unpack beta = param

    green = GreenFunc.Green2DLR{Float64}(:polarization,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1,1,ki,ni] = pifunc(k, n, param)
        end
    end
    green.dynamic=green_dyn
    return green
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    beta = 1e4
    param = Polarization.Parameter.defaultUnit(beta,1.0)
    # println("n=0")
    # println(Polarization.Polarization0_ZeroTemp(0.0, 0, param))
    # println(Polarization.Polarization0_FiniteTemp(0.0, 0, param))
    # println(Polarization.Polarization0_ZeroTemp(1e-160, 0, param))
    # println(Polarization.Polarization0_FiniteTemp(1e-160, 0, param))
    # println(Polarization.Polarization0_ZeroTemp(1e-8, 0, param))
    # println(Polarization.Polarization0_FiniteTemp(1e-8, 0, param))
    # println(Polarization.Polarization0_ZeroTemp(1.0, 0, param))
    # println(Polarization.Polarization0_FiniteTemp(1.0, 0, param))
    # println(Polarization.Polarization0_ZeroTemp(2.0, 0, param))
    # println(Polarization.Polarization0_FiniteTemp(2.0, 0, param))
    # println("n=1")
    # println(Polarization.Polarization0_ZeroTemp(0.0, 1, param))
    # println(Polarization.Polarization0_FiniteTemp(0.0, 1, param))
    # println(Polarization.Polarization0_ZeroTemp(1e-160, 1, param))
    # println(Polarization.Polarization0_FiniteTemp(1e-160, 1, param))
    # println(Polarization.Polarization0_ZeroTemp(1e-8, 1, param))
    # println(Polarization.Polarization0_FiniteTemp(1e-8, 1, param))
    # println(Polarization.Polarization0_ZeroTemp(1.0, 1, param))
    # println(Polarization.Polarization0_FiniteTemp(1.0, 1, param))
    # println(Polarization.Polarization0_ZeroTemp(2.0, 1, param))
    # println(Polarization.Polarization0_FiniteTemp(2.0, 1, param))
    Π = Polarization.Polarization0wrapped(100*param.EF,1e-4,[0.0,0.5,1.0,2.0,10.0],param)
    println(Π.dynamic)
    println(Π.dlrGrid.n)
end

