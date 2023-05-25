using MCIntegration
using ElectronGas

const para = Parameter.rydbergUnit(1.0 / 10, 4.0, 3)
const n = 0
const k = 0.1 * para.kF

function integrand(vars, config)
    n, k, para = config.userdata
    β, me, μ = para.β, para.me, para.μ
    # f(k) = 1.0 / (exp(β * (k^2 / 2 / me - μ)) + 1)
    function f(ek)
        ek = β * ek
        if ek > 0.0
            return exp(-ek) / (exp(-ek) + 1)
        else
            return 1.0 / (exp(ek) + 1)
        end
    end
    x, θ = vars[1][1], vars[2][1]
    p = x / (1 - x)
    freq = n * 2π / β * 1im
    factor = 1 / (2π)^2 * sin(θ) * p^2 / (1 - x)^2
    Ekp = (k^2 + p^2 + 2 * k * p * cos(θ)) / 2 / me - μ
    Ep = p^2 / 2 / me - μ
    w = (f(Ekp) - f(-Ep)) / (freq - Ekp - Ep) - me / p^2
    if isfinite(w) && isfinite(factor)
        return w * factor
        #return 0.0 + imag(w * factor) * 1im
    else
        return 0.0 * 1im
    end
end

function ladder(n::Int, k::Float64, para)
    return integrate(integrand;
        # var=(Continuous(2 * para.kF, 1000 * para.kF), Continuous(0.0, π * 1.0)),
        var=(Continuous(0.0, 1.0 - 1e-6), Continuous(0.0, π * 1.0)),
        dof=[[1, 1],],
        userdata=(n, k, para),
        type=ComplexF64,
        print=1,
        neval=1e6
    )
end

println(ladder(n, k, para))