using ElectronGas, JLD2
using ElectronGas.CompositeGrids
using ElectronGas.Parameters

function Ws1(q, n, param)
    if abs(q) < 1e-16
        q = 1e-16
    end
    ωn = 2 * n * π / param.β
    return 4π * param.e0^2 / (q^2 + param.Λs) * ωn^2 / (ωn^2 + param.ωp^2 / 2)
end

function Ws_freq(ωn, param)
    # ωn = 2 * n * π / param.β
    return ωn^2 / (ωn^2 + param.ωp^2 / 2)
end

function Πs0(ωn, param; ω_c=0.1param.EF)
    @unpack me, β, μ, kF, EF = param
    return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    # if abs(ωn) > ω_c
    #     return 0.0
    # else
    #     return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    # end
end

θ, rs, dim = 0.000001, 1.919, 3
param = ElectronGas.Parameter.rydbergUnit(θ, rs, dim)
wc1, wc2 = 0.1param.EF, 0.01param.EF

Ec = 100param.EF
alpha = 0.882
wgrid = CompositeGrids.CompositeG.LogDensedGrid(:gauss, [alpha / param.β, Ec], [alpha / param.β,], 24, alpha / param.β / 2, 8)
dataj = zeros(Float64, length(wgrid))
datai = zeros(Float64, length(wgrid))
for i in 1:length(wgrid)
    for j in 1:length(wgrid)
        wi, wj = wgrid[i], wgrid[j]
        dataj[j] = Πs0(wi, param; ω_c=wc1) * Πs0(wj, param; ω_c=wc2) * Ws_freq(wi - wj, param)
    end
    datai[i] = Interp.integrate1D(dataj, wgrid)
end
result = Interp.integrate1D(datai, wgrid) / π^2
println(result)



# nec = 5000
# ngrid = [n for n in -nec:nec]

# let result = 0.0
#     for i in ngrid
#         for j in ngrid
#             wi, wj = π * (2i + 1) / param.β, π * (2j + 1) / param.β
#             result += Πs0(wi, param; ω_c=wc1) * Πs0(wj, param; ω_c=wc2) * Ws_freq(wi - wj, param) / param.β^2
#         end
#     end
#     println(result)
# end