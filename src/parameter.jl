"""
Template of parameter. A submodule of ElectronGas.

    Use the convention where ħ=1, k_B=1.
    Only stores parameters that might change for purposes.
"""
module Parameter

# using Parameters
using ..Parameters

@with_kw struct Para
    WID::Int = 1

    dim::Int = 3    # dimension (D=2 or 3, doesn't work for other D!!!)
    spin::Int = 2  # number of spins

    # prime parameters
    ϵ0::Float64 = 1 / (4π)
    e0::Float64 = sqrt(2) # electron charge
    me::Float64 = 0.5  # electron mass
    EF::Float64 = 1.0     #kF^2 / (2me)
    beta::Float64 = 200
    μ::Float64 = 1.0

    # derived parameters

    # artificial parameters
    Λs::Float64 = 0.0
    Λa::Float64 = 0.0
    espin::Float64 = 0.0
end

function Base.getproperty(obj::Para, sym::Symbol)
    if sym === :β
        return obj.beta / obj.EF
    elseif sym === :n
        return (obj.dim == 3) ? (obj.EF * 2 * obj.me)^(3 / 2) / (6π^2) * obj.spin : obj.me * obj.EF / π
    elseif sym === :Rs
        return (obj.dim == 3) ? (3 / (4π * obj.n))^(1 / 3) : sqrt(1 / (πobj.n))
    elseif sym === :a0
        return 4π * obj.ϵ0 / (obj.me * obj.e0^2)
    elseif sym === :rs
        return obj.Rs / obj.a0
    elseif sym === :kF
        return sqrt(2 * obj.me * obj.EF)
    elseif sym === :e0s
        return obj.e0
    elseif sym === :e0a
        return obj.espin
    else # fallback to getfield
        return getfield(obj, sym)
    end
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
    return Para(
        dim = dim,
        spin = spin,
        ϵ0 = ϵ0,
        e0 = e0,
        me = me,
        EF = EF,
        beta = beta,
        μ = EF,
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
    ϵ0 = 1 / (4π)
    e0 = (dim == 3) ? sqrt(2 * rs / (9π / (2spin))^(1 / 3)) : sqrt(sqrt(2) * rs)
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
    ϵ0 = 1 / (4π)
    e0 = sqrt(2)
    me = 0.5
    kF = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
    EF = kF^2 / (2me)
    beta = β * EF
    return fullUnit(ϵ0, e0, me, EF, beta, dim, spin)
end

export Para, Param

end
