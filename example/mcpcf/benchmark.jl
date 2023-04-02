
# benchmark with mc results
using MCIntegration
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids
using Test

function test_rpa_interaction(t, rs)
    param = Parameter.defaultUnit(t, rs, 3)
    mint = 0.5
    minK, maxK = 0.5 * sqrt(param.T * param.me), 5param.kF
    order = 3
    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = ElectronGas.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)
    rpad, rpai = ElectronGas.Interaction.RPAwrapped(100 * param.EF, 1e-10, kgrid, param)

    tgrid = rpad.mesh[2]
    for (iq, q) in enumerate(kgrid)
        for (in, n) in enumerate(tgrid)
            n = tgrid.grid[in]
            println("q=$q, n=$n")
            W1 = Interaction.RPA(q, n, param)[1]
            println(W1)
            V = 1 / rpai[1, 1, iq]
            W2 = V * rpad[1, in, iq]
            println(W2)
            println(isapprox(W1, W2, rtol=1e-6))
        end
    end

end

function integrand_instant(p, param; k=param.kF)
    cut = 1e-16
    if abs(k - p) < cut
        return 0.0
    end
    β = param.β
    ε = p^2 / 2 / param.me - param.μ
    return p / k / (4 * π^2) * param.e0s^2 / param.ϵ0 * log(abs((k + p) / (k - p))) / (2ε) * tanh(ε * β / 2)
end

function benchmark_instant_oneloop(θ, rs, k)
    # param = Parameter.rydbergUnit(θ, rs, 3)
    param = Parameter.defaultUnit(θ, rs, 3)

    X = Continuous(0.0, 10.0)
    res = integrate((x, c) -> integrand_instant(x[1], param; k=k), solver=:vegasmc; var=X)
    report(res)
end

function integrand(idx, vars, config)
    if idx == 1
        return 1.0
    else
        norm = config.normalization
        therm = 10
        norm = sqrt(norm^2 + therm^2)

        ExtK, P = vars
        param, kgrid, data = config.userdata
        k = kgrid[ExtK[1]]
        p = P[1]

        R = Interp.linear1D(data, kgrid, p) / norm
        return -integrand_instant(p, param; k=k) * R
    end
end

function measure(idx, vars, obs, weight, config)
    config.userdata[3][vars[1][1]] += weight[1]
    if idx == 1
        obs[1][vars[1][1]] += weight[1]
    else
        obs[2][vars[1][1]] += weight[1]
    end
end

function benchmark_instant_scf(θ, rs; steps=1e6)
    param = Parameter.defaultUnit(θ, rs, 3)
    minK, maxK = 0.5 * sqrt(param.T * param.me), 10param.kF
    order = 4
    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    ExtK = Discrete(1, length(kgrid); adapt=false)
    P = Continuous(0.0, maxK; alpha=3.0, adapt=true)
    dof = [[1, 0], [1, 1]]
    data = zeros(Float64, length(kgrid))
    obs = [data, data]

    result = integrate(integrand; measure=measure, userdata=(param, kgrid, data),
        var=(ExtK, P), dof=dof, obs=obs, solver=:mcmc, neval=steps,
        print=0, block=8)
    report(result)
    println(kgrid)
    println(result.mean)
    println(result.mean[1] .+ result.mean[2])
    println(sqrt.(result.stdev[1] .^ 2 .+ result.stdev[2] .^ 2))
    # println(data)
    return result
end

test_rpa_interaction(0.2, 0.2)
# benchmark_instant_scf(0.1, 0.1; steps=1e7)

# benchmark_instant_oneloop(0.1, 3.0, 1e-8)

# function oneloop(dim, θ, rs, channel; kwargs...)
#     param = Parameter.defaultUnit(θ, rs, dim)
#     if haskey(kwargs, :int_type) && kwargs[:int_type] == :none
#         param = Parameter.Para(param; gs=0, ga=0)
#     end
#     uid = 0
#     if haskey(kwargs, :uid)
#         uid = kwargs[:uid]
#     end

#     println("dim=$dim, θ=$θ, rs=$rs, channel=$channel:")

#     lamu, R_freq, F_freq, R_ins = BSeq.linearResponse(param, channel;
#         Ntherm=0, Nmax=0, kwargs...)

#     return R_freq, F_freq, R_ins
# end

# Rw, Fw, Ri = oneloop(3, 0.1, 3.0, 0; α=0.0)

# println(Ri.mesh[2])
# println(Ri.data)