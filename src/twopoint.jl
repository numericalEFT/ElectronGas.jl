"""
Provide N-body response and correlation functions
"""
module TwoPoint

export freePropagatorT, freePropagatorΩ
export freePolarizationT
using Lehmann
using Cuba

"""
    freePropagatorT(type, τ, ω, β)

Imaginary-time propagator.

# Arguments
- `type`: symbol :fermi, :bose
- `τ`: the imaginary time, must be (-β, β]
- `ω`: dispersion ϵ_k-μ
- `β = 1.0`: the inverse temperature 
"""
@inline function freePropagatorT(type, τ, ω, β)
    return kernelT(type, τ, ω, β)
end

"""
    freePropagatorΩ(type, n, ω, β=1.0)

Matsubara-frequency kernel of different type

# Arguments
- `type`: symbol :fermi, :bose, :corr
- `n`: index of the Matsubara frequency
- `ω`: dispersion ϵ_k-μ
- `β`: the inverse temperature 
"""
@inline function freePropagatorΩ(type, n::Int, ω, β)
    return kernelΩ(type, n, ω, β)
end

@inline function freeFermiDoS(dim, kF, m, spin)
    if dim == 3
        return spin * m * kF / 2 / π^2
    else
        error("Dimension $dim not implemented!")
        # return spin/4/π
    end
end

# function LindhardΩn(dim, q, ω, β, kF, m, spin)
#     q < 0.0 && (q = -q) # Lindhard function is q symmetric

#     q2 = q^2
#     kFq = 2kF * q
#     ωn = 2π * n / β
#     D = 1 / (8kF * q)
#     NF = freeFermiDoS(dim, kF, m, spin)

#     # if ωn<=20*(q2+kFq)/(2m)
#         # careful for small q or large w
#     iw = ωn * im
#     wmq2 = iw * 2m - q^2
#     wpq2 = iw * 2m + q^2
#     C1 = log(wmq2 - kFq) - log(wmq2 + kFq)
#     C2 = log(wpq2 - kFq) - log(wpq2 + kFq)
#     res = real(-NF / 2 * (1 - D * (wmq2^2 / q^2 - 4 * kF^2) * C1 + D * (wpq2^2 / q^2 - 4 * kF^2) * C2))
#     # else
#     #     b2 = q2 * ( q2 + 12/5 * kF^2 )
#     #     c = 2*EF*kF*q2/(3*pi**2)
#     #     res = -c/(w**2 + b2)
#     # end
#     return res
# end

"""
    LindhardΩnFiniteTemperature(dim::Int, q::T, n::Int, μ::T, kF::T, β::T, m::T, spin) where {T <: AbstractFloat}

Compute the polarization function of free electrons at a given frequency. Relative Accuracy is about ~ 1e-6

# Arguments
- `dim`: dimension
- `q`: external momentum, q<1e-4 will be treated as q=0 
- `n`: externel Matsubara frequency, ωn=2π*n/β
- `μ`: chemical potential
- `kF`: Fermi momentum 
- `β`: inverse temperature
- `m`: mass
- `spin` : number of spins
"""
@inline function LindhardΩnFiniteTemperature(dim::Int, q::T, n::Int, μ::T, kF::T, β::T, m::T, spin) where {T <: AbstractFloat}
    if q < 0.0
        q = -q
    end

    if q / kF < 1.0e-10
        q = 1.0e-10 * kF
    end

    function polar(k)
        phase = T(1.0)
        if dim == 3
            phase *= k^2 / (4π^2)
        else
            error("not implemented")
        end
        ω = 2π * n / β
        ϵ = β * (k^2 - μ) / (2m)

        p = phase * Spectral.fermiDirac(ϵ) * m / k / q * log(((q^2 - 2k * q)^2 + 4m^2 * ω^2) / ((q^2 + 2k * q)^2 + 4m^2 * ω^2)) * spin

        if isnan(p)
            println("warning: integrand at ω=$ω, q=$q, k=$k is NaN!")
        end
        # println(p)
        return p
    end

    function integrand(x, f)
        # x[1]:k
        f[1] = polar(x[1] / (1 - x[1])) / (1 - x[1])^2
    end

    # TODO: use CompositeGrids to perform the integration
    result, err = Cuba.cuhre(integrand, 2, 1, rtol=1.0e-10) 
    # result, err = Cuba.vegas(integrand, 1, 1, rtol=rtol)
return result[1], err[1]
end


"""
3D Imaginary-time effective interaction derived from random phase approximation. 
Return ``dW_0(q, τ)/v_q`` where ``v_q`` is the bare Coulomb interaction and ``dW_0`` is the dynamic part of the effective interaction. 

The total effective interaction can be recoverd using,  
```math
W_0(q, τ) = v_q δ(τ) + dW_0(q, τ).
```

The dynamic contribution is the fourier transform of,
```math
dW_0(q, iω_n)=v_q^2 Π(q, iω_n)/(1-v_q Π(q, iω_n))
```
Note that this dynamic contribution ``dW_0'' diverges at small q. For this reason, this function returns ``dW_0/v_q``

# Arguments
- `vqinv`: inverse bare interaction as a function of q
- `qgrid`: one-dimensional array of the external momentum q
- `τgrid`: one-dimensional array of the imaginary-time
- `dim`: dimension
- `μ`: chemical potential
- `kF`: Fermi momentum 
- `β`: inverse temperature
- `spin` : number of spins
- `mass`: mass
"""
function dWRPA(vqinv, qgrid, τgrid, dim, μ, kF, β, spin, mass)
    @assert all(qgrid .!= 0.0)
    EF = kF^2 / (2mass)
    dlr = DLR.DLRGrid(:corr, 10EF, β, 1e-10) # effective interaction is a correlation function of the form <O(τ)O(0)>
    Nq, Nτ = length(qgrid), length(τgrid)
    Π = zeros(Complex{Float64}, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    dW0norm = similar(Π)
    for (ni, n) in enumerate(dlr.n)
        for (qi, q) in enumerate(qgrid)
            Π[qi, ni] = LindhardΩnFiniteTemperature(dim, q, n, μ, kF, β, mass, spin)[1]
        end
        dW0norm[:, ni] = @. Π[:, ni] / (vqinv - Π[:, ni])
        # println("ω_n=2π/β*$(n), Π(q=0, n=0)=$(Π[1, ni])")
        # println("$ni  $(dW0norm[2, ni])")
    end
    dW0norm = DLR.matfreq2tau(:corr, dW0norm, dlr, τgrid, axis=2) # dW0/vq in imaginary-time representation, real-valued but in complex format
    
    # println(dW0norm[1, :])
    # println(DLR.matfreq2tau(:corr, dW0norm[1, :], dlr, τgrid, axis=1))
    # println(DLR.matfreq2tau(:corr, dW0norm[2, :], dlr, τgrid, axis=1))
    # coeff = DLR.matfreq2dlr(:corr, dW0norm[1, :], dlr)
    # fitted = DLR.dlr2matfreq(:corr, coeff, dlr, dlr.n)
    # for (ni, n) in enumerate(dlr.n)
    #     println(dW0norm[1, ni], " vs ", fitted[ni])
    # end
    return real.(dW0norm) 
end

function interactionDynamic(config, qd, τIn, τOut)
    para = config.para

    dτ = abs(τOut - τIn)

    kDiQ = sqrt(dot(qd, qd))
    vd = 4π * e0^2 / (kDiQ^2 + lambda1)
        if kDiQ <= para.qgrid.grid[1]
        q = para.qgrid.grid[1] + 1.0e-6
        wd = vd * Grid.linear2D(para.dW0, para.qgrid, para.τgrid, q, dτ)
        # the current interpolation vanishes at q=0, which needs to be corrected!
    else
        wd = vd * Grid.linear2D(para.dW0, para.qgrid, para.τgrid, kDiQ, dτ) # dynamic interaction, don't forget the singular factor vq
    end

    return vd / β, wd
end

function vertexDynamic(config, qd, qe, τIn, τOut)
    vd, wd = interactionDynamic(config, qd, τIn, τOut)
    ve, we = interactionDynamic(config, qe, τIn, τOut)

    return -vd, -wd, ve, we
end

end