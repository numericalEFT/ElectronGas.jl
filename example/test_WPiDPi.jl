using ElectronGas, JLD2
using ElectronGas.CompositeGrids
using ElectronGas.Parameters

const wdratio = 0.001 # default

function Ws_freq(ωn, param)
    # ωn = 2 * n * π / param.β
    # return ωn^2 / (ωn^2 + param.ωp^2 / 2)
    return ωn^2 / (ωn^2 + 0.01 * param.EF^2)
end

function D_freq(ωn, param; ωD=wdratio * param.EF)
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

function Πr0(ωn, param; ω_c=0.1param.EF)
    @unpack me, β, μ, kF, EF = param
    if abs(ωn) < ω_c
        return 0.0
    else
        return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    end
end

function WPiDPi(θ, ωc, ωD)
    rs, dim = 1.919, 3
    param = ElectronGas.Parameter.rydbergUnit(θ, rs, dim)
    nec = floor(Int, 2param.EF / 2π * param.β + 1)
    ngrid = [n for n in -nec:nec]
    result = 0.0
    for i in ngrid
        for j in ngrid
            wi, wj = π / param.β * (2i + 1), π / param.β * (2j + 1)
            # result += Ws_freq(wi - π / param.β, param) * Πs0(wi, param; ω_c=ωc) * Πs0(wj, param; ω_c=ωc) * D_freq(wi - wj, param; ωD=ωD)
            # result += Πr0(wi, param; ω_c=ωc) * Πs0(wj, param; ω_c=ωc) * D_freq(wi - wj, param; ωD=ωD)
            result += Πs0(wj, param; ω_c=ωc) * D_freq(wi - wj, param; ωD=ωD)
        end
    end
    result = result / param.β^2
    return result
end

wdlist = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]
θ = 0.002
ωc = 0.1
WPDPs = [WPiDPi(θ, ωc, wd) for wd in wdlist]
println(WPDPs)

