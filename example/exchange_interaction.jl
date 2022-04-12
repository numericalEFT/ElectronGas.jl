"""
In this demo, we analysis the contribution of the exchange KO interaction to the Landau parameters
"""

using ElectronGas, Parameters
using Lehmann, GreenFunc, CompositeGrids

const rs = 7.0
const beta = 25.0
const mass2 = 1e-8
const dim = 3
const z = 1.0

Fp = -1.5
Fm = -0.0
# U = -0.456

U = 0.0
Cp, Cm = U / 2, -U / 2
massratio = 1.0

para = Parameter.rydbergUnit(1 / beta, rs, dim; Λs = mass2)
println(para)
const kF = para.kF
const EF = para.EF
const β = para.β
const spin = para.spin
const NF = para.NF * massratio
println(NF)

function exchange_interaction(Fp, Fm, massratio, Cp = 0.0, Cm = 0.0)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi], Wm[qi] = Interaction.KO_total(q, 0, para;
            pifunc = Polarization.Polarization0_ZeroTemp_Quasiparticle,
            landaufunc = Interaction.landauParameterConst,
            Vinv_Bare = Interaction.coulombinv,
            counter_term = Interaction.countertermConst,
            Fs = -Fp, Fa = -Fm, Cs = -Cp, Ca = -Cm, massratio = massratio)
        # instantS[qi] = Interaction.coulombinv(q, para)[1]
        # println(q, " -> ", Wp[qi] * NF, ", ", Wm[qi] * NF)
    end
    Wp *= -NF * z^2 # additional minus sign because the interaction is exchanged
    Wm *= -NF * z^2
    return Wp, Wm, θgrid
end

function Legrendre(l, func, θgrid)
    if l == 0
        return Interp.integrate1D(func .* sin.(θgrid.grid), θgrid) / 2
    elseif l == 1
        return Interp.integrate1D(func .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
    else
        error("not implemented!")
    end
end

# exchange interaction (Ws + Wa \sigma\sigma)_ex to a direct interaction Ws'+Wa' \sigma\sigma 
# # exchange S/A interaction projected to the spin-symmetric and antisymmetric parts
function exchange2direct(Wse, Wae)
    Ws = (Wse + 3 * Wae) / 2
    Wa = (Wse - Wae) / 2
    return Ws, Wa
end

function projected_exchange_interaction(l, Fp, Fm, massratio, verbose = 1)
    println("l=$l:")
    Wse, Wae, θgrid = exchange_interaction(Fp, Fm, massratio)
    Wse0 = Legrendre(l, Wse, θgrid)
    Wae0 = Legrendre(l, Wae, θgrid)
    verbose > 1 && println("Wse_l=$l=", Wse0)
    verbose > 1 && println("Wae_l=$l=", Wae0)

    Ws0, Wa0 = exchange2direct(Wse0, Wae0)
    verbose > 0 && println("Ws_l=$l=", Ws0)
    verbose > 0 && println("Wa_l=$l=", Wa0)
    return Ws0, Wa0
end

Ws0, Wa0 = projected_exchange_interaction(0, Fp, Fm, massratio)
Wsc, Wac = exchange2direct(Fp, Fm)
println(Ws0 + Wsc)
println(Wa0 + Wac)
# projected_exchange_interaction(1, Fp, Fm, massratio)
exit(0)

Fs, Fa = Fp, Fm
mix = 0.2
for i = 1:100
    println("iteration: $i")
    global Fs, Fa
    nFs, nFa = projected_exchange_interaction(0, Fs, Fa, massratio)
    Fs = Fs * (1 - mix) + 2 * nFs * mix
    Fa = Fa * (1 - mix) + nFa * mix
    Fa = 0.0
end
println("Fs = ", Fs)
println("Fa = ", Fa)

# Ws1 = Interp.integrate1D(Ws .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
# Wa1 = Interp.integrate1D(Wa .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2

# println("l=1:")
# println("F1+=", Ws1)
# println("F1-=", Wa1)
# println("m^*/m = ", 1 + Ws1)

# println(sqrt(4π * para.e0^2 * para.NF) / para.kF)
