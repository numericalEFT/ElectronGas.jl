using ElectronGas
using JLD2

function load_AB()
end

function init_AB()
end

function Πs0(ωn, param; ω_c=0.02param.EF)
    @unpack me, β, μ, kF, EF = param
    return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
end
