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
    β::Float64 = 200 # bare inverse temperature
    μ::Float64 = 1.0

    # artificial parameters
    Λs::Float64 = 0.0   # Yukawa-type spin-symmetric interaction  ~1/(q^2+Λs)
    Λa::Float64 = 0.0   # Yukawa-type spin-asymmetric interaction ~1/(q^2+Λa)
    espin::Float64 = 0.0
end

function Base.getproperty(obj::Para, sym::Symbol)
    if sym === :beta # dimensionless beta
        return obj.β * obj.EF
    elseif sym === :Θ # dimensionless temperature
        return 1.0 / obj.β / obj.EF
    elseif sym === :T
        return 1.0 / obj.β
    elseif sym === :n
        return (obj.dim == 3) ? (obj.EF * 2 * obj.me)^(3 / 2) / (6π^2) * obj.spin : obj.me * obj.EF / π
    elseif sym === :Rs
        return (obj.dim == 3) ? (3 / (4π * obj.n))^(1 / 3) : sqrt(1 / (π * obj.n))
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
    function fullUnit(ϵ0, e0, me, EF, β)

generate Para with a complete set of parameters, no value presumed.

#Arguments:
 - ϵ0: vacuum permittivity
 - e0: electron charge
 - me: electron mass
 - EF: Fermi energy
 - β: inverse temperature
"""
@inline function fullUnit(ϵ0, e0, me, EF, β, dim = 3, spin = 2; kwargs...)
    para = Para(dim = dim,
        spin = spin,
        ϵ0 = ϵ0,
        e0 = e0,
        me = me,
        EF = EF,
        β = β,
        μ = EF
    )
    return reconstruct(para, kwargs...)
end

"""
    function defaultUnit(Θ, rs)

assume 4πϵ0=1, me=0.5, EF=1

#Arguments:
 - Θ: dimensionless temperature. Since EF=1 we have β=beta
 - rs: Wigner-Seitz radius over Bohr radius.
"""
@inline function defaultUnit(Θ, rs, dim = 3, spin = 2; kwargs...)
    ϵ0 = 1 / (4π)
    e0 = (dim == 3) ? sqrt(2 * rs / (9π / (2spin))^(1 / 3)) : sqrt(sqrt(2) * rs)
    me = 0.5
    EF = 1
    β = 1 / Θ / EF
    return fullUnit(ϵ0, e0, me, EF, β, dim, spin; kwargs...)
end


"""
    function rydbergUnit(Θ, rs)

assume 4πϵ0=1, me=0.5, e0=sqrt(2)

#Arguments:
 - Θ: dimensionless temperature. beta could be different from β
 - rs: Wigner-Seitz radius over Bohr radius.
"""
@inline function rydbergUnit(Θ, rs, dim = 3, spin = 2; kwargs...)
    ϵ0 = 1 / (4π)
    e0 = sqrt(2)
    me = 0.5
    kF = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
    EF = kF^2 / (2me)
    β = 1 / Θ / EF
    return fullUnit(ϵ0, e0, me, EF, β, dim, spin; kwargs...)
end

export Para, Param

end
