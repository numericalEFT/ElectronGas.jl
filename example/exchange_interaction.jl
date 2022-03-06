"""
In this demo, we analysis the contribution of the exchange KO interaction to the Landau parameters
"""

using ElectronGas, Parameters
using Lehmann, GreenFunc, CompositeGrids

const rs = 5.0
const beta = 1000.0
const mass2 = 1e-8
const dim = 3

Fp = -0.0
Fm = -0.0
massratio = 0.95

para = Parameter.rydbergUnit(1 / beta, rs, dim; Λs = mass2)
println(para)
const kF = para.kF
const EF = para.EF
const β = para.β
const spin = para.spin
const NF = para.NF * massratio
println(NF)

θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

Wp = zeros(Float64, length(qs))
Wm = zeros(Float64, length(qs))
for (qi, q) in enumerate(qs)
    Wp[qi], Wm[qi] = Interaction.KO_total(q, 0, para;
        pifunc = Polarization.Polarization0_ZeroTemp_Quasiparticle,
        landaufunc = Interaction.landauParameterConst,
        Vinv_Bare = Interaction.coulombinv,
        counter_term = Interaction.counterterm,
        Fs = -Fp, Fa = -Fm, massratio = massratio)
    # instantS[qi] = Interaction.coulombinv(q, para)[1]
    # println(q, " -> ", Ws[qi] * NF, ", ", Wa[qi] * NF)
end
Wp *= NF
Wm *= NF

# exchange KO interaction projected to the spin-symmetric and antisymmetric parts
Ws, Wa = -(Wp + 3 * Wm) / 2, -(Wp - Wm) / 2

Ws0 = Interp.integrate1D(Ws .* sin.(θgrid.grid), θgrid) / 2
Wa0 = Interp.integrate1D(Wa .* sin.(θgrid.grid), θgrid) / 2
println("l=0:")
println("F0+=", Ws0)
println("F0-=", Wa0)

Ws1 = Interp.integrate1D(Ws .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
Wa1 = Interp.integrate1D(Wa .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2

println("l=1:")
println("F1+=", Ws1)
println("F1-=", Wa1)
println("m^*/m = ", 1 + Ws1)

println(sqrt(4π * para.e0^2 * para.NF) / para.kF)
