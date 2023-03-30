using MCIntegration
using LegendrePolynomials

const steps = 1e6
const ℓ = 0


include("propagators.jl")
using .Propagators
using .Propagators: G0, interaction, response

const param = Propagators.Parameter.defaultUnit(0.1, 2.0)

function integrand(vars, config)
    norm = config.normalization
    if norm == 0.0
        norm = 1.0
    end

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
    G1 = G0(t1, p, funcs)
    G021 = G0(t1, -p, funcs)
    G022 = G0(t2, -p, funcs)
    R0 = response(p, funcs; norm=norm)
    R = response(t1 - t2, p, funcs; norm=norm)

    result1 = -p^2 / (4π^2) * PLX * V * G1 * (G021 * R0 + G022 * R)

    W = interaction(t, q, funcs)
    G21 = G0(t1 - t, -p, funcs)
    G22 = G0(t2 - t, -p, funcs)

    result2 = -p^2 / (4π^2) * PLX * W * G1 * (G21 * R0 + G22 * R)
    if isnan(result1) || isnan(result2) || isinf(result1) || isinf(result2)
        println("t=$t, k=$k, x=$x, t1=$t1, t2=$t2, p=$p")
        println("PLX=$PLX, q=$q, V=$V,G1=$G1, G021=$G021, G022=$G022, R0=$R0, R=$R")
        println("W=$W, G21=$G21, G22=$G22")
        error("nan!")
    end
    return result1 / length(extT), result2
end

function measure(vars, obs, weight, config)
    extt, extk = vars[1], vars[2]
    # obs[1][extk[1]] += weight[1]
    # obs[2][extt[1], extk[1]] += weight[2]
    funcs = config.userdata
    Ri, Rt = funcs.Ri, funcs.Rt
    Ri.data[extk[1]] += weight[1]
    Rt.data[extt[1], extk[1]] += weight[2]
    obs[1] = config.normalization^2
end

function run(steps, param, alg=:vegas)
    println("Prepare propagators")
    mint = 0.01
    minK, maxK = 0.01 * sqrt(param.T * param.me), 10param.kF
    order = 4

    rpai, rpat = Propagators.rpa(param; mint=mint, minK=minK, maxK=maxK, order=order)
    Ri, Rt = Propagators.initR(param; mint=mint, minK=minK, maxK=maxK, order=order)
    println(size(Rt))

    # userdata
    funcs = Propagators.Funcs(param, rpai, rpat, Ri, Rt)

    println("Prepare variables")
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    T = Continuous(0.0, param.β; alpha=3.0, adapt=true)
    P = Continuous(extK.bound[1], extK.bound[2]; alpha=3.0, adapt=true)
    X = Continuous(-1.0, 1.0; alpha=3.0, adapt=true)

    ExtT = Discrete(1, length(extT); adapt=false)
    ExtK = Discrete(1, length(extK); adapt=false)

    # ExtT, ExtK, X, T, P
    dof = [[0, 1, 1, 2, 1], [1, 1, 1, 2, 1]]
    obs = [0.0, 0.0]

    println("Start")
    result = integrate(integrand; measure=measure, userdata=funcs,
        var=(ExtT, ExtK, X, T, P), dof=dof, obs=obs, solver=alg,
        neval=steps, print=0, block=8, type=ComplexF64)

    norm = result[1][1]
    funcs.Ri.data .= funcs.Ri.data ./ norm
    funcs.Rt.data .= funcs.Rt.data ./ norm

    return result, funcs
end

if abspath(PROGRAM_FILE) == @__FILE__
    result, funcs = run(steps, param, :vegasmc)
    Ri, Rt = funcs.Ri, funcs.Rt
    println(Ri.mesh[1])
    println(Ri.data)
    println("R0=$(Propagators.R0(Ri, Rt, param))")
end