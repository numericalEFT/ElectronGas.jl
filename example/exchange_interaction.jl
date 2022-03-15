"""
In this demo, we analysis the contribution of the exchange KO interaction to the Landau parameters
"""

using ElectronGas, Parameters
using Lehmann, GreenFunc, CompositeGrids

const rs = 5.3
const beta = 1000.0
const mass2 = 1e-8
const dim = 3

Fp = -1.0
Fm = -0.5
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
    Wp *= -NF # additional minus sign because the interaction is exchanged
    Wm *= -NF
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
function exchange2direct(Wse, Wae)
    Ws = (Wse + 3 * Wae) / 2
    Wa = (Wse - Wae) / 2
    return Ws, Wa
end

# Wp0 = Interp.integrate1D(Wp .* sin.(θgrid.grid), θgrid) / 2
# Wm0 = Interp.integrate1D(Wm .* sin.(θgrid.grid), θgrid) / 2
println("l=0:")
# println("F0+=", Ws0 + Fp)
# println("F0-=", Wa0 + Fm)
Wse, Wae, θgrid = exchange_interaction(Fp, Fm, massratio)
Wse0 = Legrendre(0, Wse, θgrid)
Wae0 = Legrendre(0, Wae, θgrid)
println("Wse_l=0=", Wse0)
println("Wae_l=0=", Wae0)

Ws0, Wa0 = exchange2direct(Wse0, Wae0)
println("Ws_l=0=", Ws0)
println("Wa_l=0=", Wa0)

# # exchange KO interaction projected to the spin-symmetric and antisymmetric parts
# Ws, Wa = -(Wp + 3 * Wm) / 2, -(Wp - Wm) / 2

# Ws0 = Interp.integrate1D(Ws .* sin.(θgrid.grid), θgrid) / 2
# Wa0 = Interp.integrate1D(Wa .* sin.(θgrid.grid), θgrid) / 2
# println("l=0:")
# # println("F0+=", Ws0 + Fp)
# # println("F0-=", Wa0 + Fm)
# println("F0+ - Fp=", Ws0)
# println("F0- - Fm=", Wa0)

# Ws1 = Interp.integrate1D(Ws .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
# Wa1 = Interp.integrate1D(Wa .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2

# println("l=1:")
# println("F1+=", Ws1)
# println("F1-=", Wa1)
# println("m^*/m = ", 1 + Ws1)

# println(sqrt(4π * para.e0^2 * para.NF) / para.kF)
