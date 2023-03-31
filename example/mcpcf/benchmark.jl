
# benchmark with mc results
using MCIntegration
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

param = Parameter.rydbergUnit(0.1, 3.0, 3)
# param = Parameter.defaultUnit(0.1, 3.0, 3)

function integrand(p; param=param, k=param.kF)
    cut = 1e-16
    if abs(k - p) < cut
        return 0.0
    end
    β = param.β
    ε = p^2 / 2 / param.me - param.μ
    return p / k / (4 * π^2) * param.e0s^2 / param.ϵ0 * log(abs((k + p) / (k - p))) / (2ε) * tanh(ε * β / 2)
end

X = Continuous(0.0, 10.0)
res = integrate((x, c) -> integrand(x[1]; k=1e-8), solver=:vegasmc; var=X)
# res = integrate((x, c) -> integrand(x[1]; k=1), solver=:vegasmc; var=X)
# res = integrate((x, c) -> integrand(x[1]; k=10), solver=:vegasmc; var=X)
report(res)

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