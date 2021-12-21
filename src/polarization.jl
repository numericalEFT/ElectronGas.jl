module Polarization

using Parameters, GreenFunc, CompositeGrids

export Polarization0_ZeroTemp

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

include(srcdir*"/parameter.jl")
using .Parameter

include(srcdir*"/convention.jl")
using .Convention

function ΠT_integrand(k, q, ω, param)
    @unpack me, beta, μ, kF = param
     return k*me/(2*π^2*q)/(exp(beta*(k^2/2/me-μ))+1)*log1p((8*k*q^3)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    # if ω==0 && abs(q*(k*2-q)*beta)<1e-16
    #     return 0.0
    # elseif ω!=0 && k^2*q^2/4/me^2<ω^2*1e-9
    #     2*q^2*k^2/(π^2*ω^2)/(exp(beta*(k^2/2/me-μ))+1)/(8*me)
    # else
        # return k*me/(2*π^2*q)/(exp(beta*(k^2/2/me-μ))+1)*log((4*me^2*ω^2+(q^2+2*k*q)^2)/(4*me^2*ω^2+(q^2-2*k*q)^2))
        # return k*me/(2*π^2*q)/(exp(beta*(k^2/2/me-μ))+1)*log1p((8*k*q^3)/(4*me^2*ω^2+(q^2-2*k*q)^2))
    # end
end

function Polarization0_FiniteTemp(q, n, param)
    @unpack me, kF, beta = param
    kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 20*kF], [0.5*q, kF], 20, 1e-6*min(q,kF), 10)
    integrand = zeros(Float64, kgrid.size)
    for (ki, k) in enumerate(kgrid.grid)
        integrand[ki] = ΠT_integrand(k, q, 2π*n/beta, param)
    end

    return Interp.integrate1D(integrand, kgrid)
end

"""
    function Polarization0_ZeroTemp(q, n, param)

Zero temperature Π0 function for matsubara frequency and momentum. For low temperature the finite temperature
polarization could be approximated with this function.
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F)

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function Polarization0_ZeroTemp(q, n, param)
    
    @unpack me, kF, rs, e0, beta , mass2, ϵ0 = param
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
    function Polarization0_ZeroTemp(T, tgrid, sgrid, param)

Zero temperature Π0 function for matsubara frequency and momentum. For low temperature the finite temperature
polarization could be approximated with this function.
Assume G_0^{-1} = iω_n - (k^2/(2m) - E_F).
Return full polarization0 function stored in GreenFunc.GreenBasic.Green2.

#Arguments:
 - T: type of data stored
 - tgrid: matsubara frequency grid
 - sgrid: momentum grid
 - param: other system parameters
"""
function Polarization0_ZeroTemp(T, tgrid::TGT, sgrid::SGT, param) where{TGT, SGT}
    green = GreenFunc.GreenBasic.Green2{T}(:freq, :mom, :bose, [0.0,], tgrid, sgrid)
    for (k, ki) in enumerate(sgrid)
        for (n, ni) in enumerate(tgrid)
            green.value[ni, ki, 1] = Polarization0_ZeroTemp(k, n, param)
        end
    end

    return green
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    beta = 1e8
    param = Polarization.Parameter.defaultUnit(beta,1.0)
    println(Polarization.Polarization0_ZeroTemp(1e-8, 0, param))
    println(Polarization.Polarization0_FiniteTemp(1e-8, 0, param))
    println(Polarization.Polarization0_ZeroTemp(1.0, 0, param))
    println(Polarization.Polarization0_FiniteTemp(1.0, 0, param))
    println(Polarization.Polarization0_ZeroTemp(2.0, 0, param))
    println(Polarization.Polarization0_FiniteTemp(2.0, 0, param))
    println(Polarization.Polarization0_ZeroTemp(1e-8, 1, param))
    println(Polarization.Polarization0_FiniteTemp(1e-8, 1, param))
    println(Polarization.Polarization0_ZeroTemp(1.0, 1, param))
    println(Polarization.Polarization0_FiniteTemp(1.0, 1, param))
    println(Polarization.Polarization0_ZeroTemp(2.0, 1, param))
    println(Polarization.Polarization0_FiniteTemp(2.0, 1, param))
end

