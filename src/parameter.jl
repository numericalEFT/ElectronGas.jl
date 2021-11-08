"""
Template of parameter.

    Use the convention where ħ=1, k_B=1.
    Only stores parameters that might change for purposes.
"""
module Parameter

using Parameters

@with_kw struct Para
    WID::Int = 1

    dim::Int = 3    # dimension (D=2 or 3, doesn't work for other D!!!)
    spin::Int = 2  # number of spins

    # prime parameters
    ϵ0::Float64 = 1/(4π)
    e0::Float64 = sqrt(2) # electron charge
    me::Float64 = 0.5  # electron mass
    EF::Float64 = 1.0     #kF^2 / (2me)
    beta::Float64 = 200

    # derived parameters
    β::Float64 = beta/EF
    n::Float64 = (dim == 3) ? (EF*2*me)^(3/2)/(6π^2)*spin : me*EF/π
    Rs::Float64 = (dim == 3) ? (3 / (4π*n))^(1 / 3) : sqrt(1/(πn))
    a0::Float64 = 4π*ϵ0/(me*e0^2)
    rs::Float64 = Rs/a0
    kF::Float64 = sqrt(2*me*EF)

    # artificial parameters
    mass2::Float64 = 0.000
end

"""
    function fullUnit(ϵ0, e0, me, EF, beta)

generate Para with a complete set of parameters, no value presumed.

#Arguments:
 - ϵ0: vacuum permittivity
 - e0: electron charge
 - me: electron mass
 - EF: Fermi energy
 - beta: inverse temperature
"""
@inline function fullUnit(ϵ0, e0, me, EF, beta, dim = 3, spin = 2)
    β::Float64 = beta/EF
    n::Float64 = (dim == 3) ? (EF*2*me)^(3/2)/(6π^2)*spin : me*EF/π
    Rs::Float64 = (dim == 3) ? (3 / (4π*n))^(1 / 3) : sqrt(1/(πn))
    a0::Float64 = 4π*ϵ0/(me*e0^2)
    rs::Float64 = Rs/a0
    kF::Float64 = sqrt(2*me*EF)
    return Para(
        dim=dim,
        spin=spin,
        ϵ0=ϵ0,
        e0=e0,
        me=me,
        EF=EF,
        beta=beta,
        β=β,
        n=n,
        Rs=Rs,
        a0=a0,
        rs=rs,
        kF=kF
    )
end

"""
    function defaultUnit(β, rs)

assume 4πϵ0=1, me=0.5, EF=1

#Arguments:
 - β: inverse temperature. Since EF=1 we have beta=β
 - rs: Wigner-Seitz radius over Bohr radius.
"""
@inline function defaultUnit(β, rs, dim = 3, spin = 2)
    ϵ0 = 1/(4π)
    e0 = (dim == 3) ? sqrt(2*rs / (9π / (2spin))^(1 / 3) ) : sqrt(sqrt(2)*rs)
    me = 0.5
    EF = 1
    beta = β * EF
    return fullUnit(ϵ0, e0, me, EF, beta, dim, spin)
end


"""
    function rydbergUnit(β, rs)

assume 4πϵ0=1, me=0.5, e0=sqrt(2)

#Arguments:
 - β: inverse temperature. Since EF=1 we have beta=β
 - rs: Wigner-Seitz radius over Bohr radius.
"""
@inline function rydbergUnit(β, rs, dim = 3, spin = 2)
    ϵ0 = 1/(4π)
    e0 = sqrt(2)
    me = 0.5
    kF = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
    EF = kF^2/(2me)
    beta = β * EF
    return fullUnit(ϵ0, e0, me, EF, beta, dim, spin)
end

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

if !@isdefined Param
    try
        include(rundir*"/para.jl")
    catch ee
        if isa(ee, LoadError)
            println("Load failed. Generating default parameters instead.")
            Param = defaultUnit(100, 1)
        end
    end
end

export Para, Param

end
