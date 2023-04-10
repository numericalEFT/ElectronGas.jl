using MCIntegration
using LegendrePolynomials

include("propagators.jl")
using .Propagators
using .Propagators: G0, interaction, response

const fname = "run/data/PCFdata_3003.jld2"
const steps = 1e7 # 2e8/hr
const ℓ = 0
const θ, rs = 0.1, 0.3
const param = Propagators.Parameter.defaultUnit(θ, rs, 3)
const α = 0.8
println(param)

function integrand(vars, config)
    norm = config.normalization
    therm = 100000
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
    # R0 = response(p, funcs; norm=norm) / param.β
    # R = response(t1 - t2, p, funcs; norm=norm)
    R0 = response(p, funcs) / param.β
    R = response(t1 - t2, p, funcs)
    # R0 = 1.0 / param.β
    # R = 0.0

    result1 = -p^2 / (4π^2) * PLX * V * G1 * (G021 * R0 + G022 * R)
    # result1 = result1 / length(extT)

    W = interaction(t, q, funcs) * V
    # W = 0.0
    G21 = G0(t1 - t, -p, funcs)
    G22 = G0(t2 - t, -p, funcs)

    result2 = -p^2 / (4π^2) * PLX * W * G1 * (G21 * R0 + G22 * R)
    # if isnan(result1) || isnan(result2) || isinf(result1) || isinf(result2)
    #     println("t=$t, k=$k, x=$x, t1=$t1, t2=$t2, p=$p")
    #     println("PLX=$PLX, q=$q, V=$V,G1=$G1, G021=$G021, G022=$G022, R0=$R0, R=$R")
    #     println("W=$W, G21=$G21, G22=$G22")
    #     error("nan!")
    # end
    return 1.0, result1, result2
end

function measure(vars, obs, weight, config)
    extt, extk = vars[1], vars[2]
    # funcs = config.userdata
    # Ri, Rt = funcs.Ri, funcs.Rt
    # Ri.data[extk[1]] += weight[1] + weight[2]
    # Rt.data[extt[1], extk[1]] += weight[3]
    obs[1][extk[1]] += weight[1]
    obs[2][extk[1]] += weight[2]
    obs[3][extt[1], extk[1]] += weight[3]
end

function run(steps, param, alg=:vegas)
    println("Prepare propagators")

    mint = 0.001
    minK, maxK = 0.001 * sqrt(param.T * param.me), 10param.kF
    order = 6
    rpai, rpat = Propagators.rpa(param; mint=mint, minK=minK, maxK=maxK, order=order)

    mint = 0.1
    minK, maxK = 0.1 * sqrt(param.T * param.me), 10param.kF
    order = 3
    Ri, Rt = Propagators.loadR(fname, param; mint=mint, minK=minK, maxK=maxK, order=order)
    # Ri, Rt = Propagators.initR(param; mint=mint, minK=minK, maxK=maxK, order=order)
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
        neval=steps, print=-1, block=8, type=ComplexF64)
    # println(result.mean)
    funcs.Ri.data .= result.mean[1] .+ result.mean[2]
    funcs.Rt.data .= result.mean[3]

    return result, funcs
end

if abspath(PROGRAM_FILE) == @__FILE__
    result, funcs = run(steps, param, :vegasmc)
    Ri, Rt = funcs.Ri, funcs.Rt
    println(Ri.mesh[1])
    println(real(Ri.data))
    # println(result[1][1])
    println("R0=$(real(Propagators.R0(Ri, Rt, param)))")
end