using MCIntegration
using LegendrePolynomials

include("propagators.jl")
using .Propagators
using .Propagators: G0, interaction, response

const steps = 1e7
const ℓ = 0
const param = Propagators.Parameter.defaultUnit(0.1, 0.1, 3)
println(param)

function integrand(idx, vars, config)
    if idx == 1
        return 1.0
    elseif idx == 2
        return integrand2(vars, config)
    elseif idx == 3
        return integrand3(vars, config)
    end
end

function integrand2(vars, config)
    norm = config.normalization
    therm = 10
    norm = sqrt(norm^2 + therm^2)

    ExtT, ExtK, X, T, P = vars
    funcs = config.userdata
    Rt = funcs.Rt
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    param = funcs.param

    t = extT[ExtT[1]]
    k = extK[ExtK[1]]
    x = X[1]
    t1, t2 = T[1], T[2]
    p = P[1]

    PLX = Pl(x, ℓ)
    q = sqrt(k^2 + p^2 + 2 * k * p * x)
    V = 1.0 / interaction(q, funcs)
    # V = coulomb(q, funcs)[1]
    G1 = G0(t1, p, funcs)
    G021 = G0(t1, -p, funcs)
    G022 = G0(t2, -p, funcs)
    R0 = response(p, funcs; norm=norm) / param.β
    R = response(t1 - t2, p, funcs; norm=norm)
    # R0 = 1.0 / param.β
    # R = 0.0

    result1 = -p^2 / (4π^2) * PLX * V * G1 * (G021 * R0 + G022 * R)
    return result1
end

function integrand3(vars, config)
    norm = config.normalization
    therm = 10
    norm = sqrt(norm^2 + therm^2)

    ExtT, ExtK, X, T, P = vars
    funcs = config.userdata
    Rt = funcs.Rt
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    param = funcs.param

    t = extT[ExtT[1]]
    k = extK[ExtK[1]]
    x = X[1]
    t1, t2 = T[1], T[2]
    p = P[1]

    PLX = Pl(x, ℓ)
    q = sqrt(k^2 + p^2 + 2 * k * p * x)
    V = 1.0 / interaction(q, funcs)
    # V = coulomb(q, funcs)[1]
    G1 = G0(t1, p, funcs)
    R0 = response(p, funcs; norm=norm) / param.β
    R = response(t1 - t2, p, funcs; norm=norm)
    # R0 = 1.0 / param.β
    # R = 0.0
    W = interaction(t, q, funcs) * V
    # W = 0.0
    G21 = G0(t1 - t, -p, funcs)
    G22 = G0(t2 - t, -p, funcs)

    result2 = -p^2 / (4π^2) * PLX * W * G1 * (G21 * R0 + G22 * R)
    return result2
end

function measure(idx, vars, obs, weight, config)
    extt, extk = vars[1], vars[2]
    funcs = config.userdata
    Ri, Rt = funcs.Ri, funcs.Rt
    if idx == 1
        obs[1][extk[1]] += weight[1]
        Ri.data[extk[1]] += weight[1]
    elseif idx == 2
        obs[2][extk[1]] += weight[1]
        Ri.data[extk[1]] += weight[1]
    elseif idx == 3
        obs[3][extt[1], extk[1]] += weight[1]
        Rt.data[extt[1], extk[1]] += weight[1]
    end
end

function run(steps, param)
    alg = :mcmc
    println("Prepare propagators")
    mint = 0.001
    minK, maxK = 0.001 * sqrt(param.T * param.me), 10param.kF
    order = 6
    rpai, rpat = Propagators.rpa(param; mint=mint, minK=minK, maxK=maxK, order=order)

    mint = 0.5
    minK, maxK = 0.5 * sqrt(param.T * param.me), 10param.kF
    order = 4
    Ri, Rt = Propagators.initR(param; mint=mint, minK=minK, maxK=maxK, order=order)
    println(size(Rt))

    # userdata
    funcs = Propagators.Funcs(param, rpai, rpat, Ri, Rt)

    println("Prepare variables")
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    T = Continuous(0.0, param.β; alpha=3.0, adapt=true)
    P = Continuous(0.0, maxK; alpha=3.0, adapt=true)
    X = Continuous(-1.0, 1.0; alpha=3.0, adapt=true)

    ExtT = Discrete(1, length(extT); adapt=false)
    ExtK = Discrete(1, length(extK); adapt=false)

    # ExtT, ExtK, X, T, P
    dof = [[0, 1, 0, 0, 0], [0, 1, 1, 2, 1], [1, 1, 1, 2, 1]]
    obs = [zeros(ComplexF64, size(Ri)), zeros(ComplexF64, size(Ri)), zeros(ComplexF64, size(Rt))]

    println("Start")
    result = integrate(integrand; measure=measure, userdata=funcs,
        var=(ExtT, ExtK, X, T, P), dof=dof, obs=obs, solver=alg,
        neval=steps, print=0, block=8, type=ComplexF64)

    println(result.mean)
    funcs.Ri.data .= result.mean[1] .+ result.mean[2]
    funcs.Rt.data .= result.mean[3]

    return result, funcs
end

if abspath(PROGRAM_FILE) == @__FILE__
    result, funcs = run(steps, param)
    Ri, Rt = funcs.Ri, funcs.Rt
    println(Ri.mesh[1])
    println(real(Ri.data))
    # println(result[1][1])
    println("R0=$(real(Propagators.R0(Ri, Rt, param)))")
end