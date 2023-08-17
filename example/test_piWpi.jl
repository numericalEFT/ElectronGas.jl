using ElectronGas, JLD2
using ElectronGas.CompositeGrids
using ElectronGas.Parameters

const wdratio = 0.001
function Ws1(q, n, param)
    if abs(q) < 1e-16
        q = 1e-16
    end
    ωn = 2 * n * π / param.β
    return 4π * param.e0^2 / (q^2 + param.Λs) * ωn^2 / (ωn^2 + param.ωp^2 / 2)
end

function Ws_freq(ωn, param)
    # ωn = 2 * n * π / param.β
    # return ωn^2 / (ωn^2 + param.ωp^2 / 2)
    return ωn^2 / (ωn^2 + 0.01 * param.EF^2)
end

function D_freq(ωn, param)
    ωD = wdratio * param.EF
    # ωn = 2 * n * π / param.β
    return ωD^2 / (ωn^2 + ωD^2)
end

function Πs0(ωn, param; ω_c=0.1param.EF)
    @unpack me, β, μ, kF, EF = param
    # return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    if abs(ωn) > ω_c
        return 0.0
    else
        return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    end
end

function calc_piwpi(θ, wc1_ratio, wc2_ratio)
    rs, dim = 1.919, 3
    param = ElectronGas.Parameter.rydbergUnit(θ, rs, dim)
    wc1, wc2 = wc1_ratio * param.EF, wc2_ratio * param.EF

    Ec = 10param.EF
    alpha = 0.882
    wgrid = CompositeGrids.CompositeG.LogDensedGrid(:gauss, [alpha / param.β, Ec], [alpha / param.β, 0.1param.EF], 30, alpha / param.β / 2, 8)
    dataj = zeros(Float64, length(wgrid))
    datai = zeros(Float64, length(wgrid))
    for i in 1:length(wgrid)
        for j in 1:length(wgrid)
            wi, wj = wgrid[i], wgrid[j]
            # wi, wj = wgrid[1], wgrid[j]
            # dataj[j] = Πs0(wi, param; ω_c=wc1) * Πs0(wj, param; ω_c=wc2) * Ws_freq(wi - wj, param)
            # dataj[j] = Ws_freq(wi - wgrid[1], param) * Ws_freq(wj - wgrid[1], param) * Πs0(wi, param; ω_c=wc2) * D_freq(wi - wj, param) * Πs0(wj, param; ω_c=wc2)
            # dataj[j] = Ws_freq(wj - wgrid[1], param) * Πs0(wj, param; ω_c=wc2)
            dataj[j] = Πs0(wi, param; ω_c=wc1) * Πs0(wj, param; ω_c=wc2) * D_freq(wi - wj, param) * Ws_freq(wi - wgrid[1], param)
        end
        datai[i] = Interp.integrate1D(dataj, wgrid)
    end
    result = Interp.integrate1D(datai, wgrid) / π^2
    return result
end

# println(calc_piwpi(1e-6, 0.1, wdratio))
num = 9
# θ1 = wdratio
θ1 = 0.001
# betas = [400 * 2^(i - 1) for i in 1:num]
betas = [1 / θ1 * 2^(i - 1) for i in 1:num]
lnbetas = [log10(beta) for beta in betas]
# piwpis = [calc_piwpi(1 / beta, 0.05, 0.05) for beta in betas]
piwpis = [calc_piwpi(1 / beta, 0.1, 0.1) for beta in betas]
# piwpis = [calc_piwpi(1 / beta, 10 / betas[1], 0.1) for beta in betas]
println(lnbetas)
println(piwpis)

println((piwpis[end] - piwpis[1]) / (lnbetas[end] - lnbetas[1]))

# using Plots
# plt = plot(lnbetas, piwpis)
# display(plt)
# readline()

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