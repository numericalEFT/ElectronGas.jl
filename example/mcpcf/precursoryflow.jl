using MCIntegration
using LegendrePolynomials

const steps = 1e5
const ℓ = 0


include("propagators.jl")
using .Propagators
using .Propagators: G0, interaction, response

const param = Propagators.Parameter.defaultUnit(0.01, 2.0)

function integrand(idx, vars, config)
    if idx == 1
        return integrand1(vars, config)
    elseif idx == 2
        return integrand2(vars, config)
    else
        error("idx=$idx invalid!")
    end
end

function integrand1(vars, config)
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
    V = interaction(q, funcs)
    G1 = G0(t1, p)
    G21 = G0(t1, -p)
    G22 = G0(t2, -p)
    R0 = response(p, funcs)
    R = response(t1 - t2, p, funcs)

    return p^2 / (4π^2) * PLX * V * G1 * (G21 * R0 + G22 * R)
end

function integrand2(vars, config)
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
    W = interaction(t, q, funcs)
    G1 = G0(t1, p)
    G21 = G0(t1 - t, -p)
    G22 = G0(t2 - t, -p)
    R0 = response(p, funcs)
    R = response(t1 - t2, p, funcs)

    return p^2 / (4π^2) * PLX * W * G1 * (G21 * R0 + G22 * R)
end

function measure(idx, vars, obs, weight, config)
    if idx == 1
        extk = vars[2]
        obs[1][extk[1]] += weight[1]
    elseif idx == 2
        extt, extk = vars[1], vars[2]
        obs[2][extt[1], extk[1]] += weight[1]
    end
end

function run(steps, param, alg=:vegas)
    rpai, rpat = Propagators.rpa(param)
    Ri, Rt = Propagators.initR(param)

    # userdata
    funcs = Propagators.Funcs(param, rpai, rpat, Ri, Rt)

    extT, extK = Rt.mesh[1], Rt.mesh[2]
    T = Continuous(0.0, param.β; alpha=3.0, adapt=true)
    P = Continuous(extK.bound[1], extK.bound[2]; alpha=3.0, adapt=true)
    X = Continuous(-1.0, 1.0; alpha=3.0, adapt=true)

    ExtT = Discrete(1, length(extT); adapt=false)
    ExtK = Discrete(1, length(extK); adapt=false)

    # ExtT, ExtK, X, T, P
    dof = [[0, 1, 1, 2, 1], [1, 1, 1, 2, 1]]
    obs = [Ri.data, Rt.data]

    result = integrate(integrand; measure=measure, userdata=funcs,
        var=(ExtT, ExtK, X, T, P), dof=dof, obs=obs, solver=alg,
        neval=steps, print=-1, block=8)

    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run(steps, param)
end