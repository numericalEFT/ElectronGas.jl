"""
In this demo, we analysis the contribution of the exchange KO interaction to the Landau parameters
"""

using ElectronGas, Parameters
using Lehmann, GreenFunc, CompositeGrids

const rs = 5.0
const beta = 1000.0
const mass2 = 1e-8
const dim = 3
θ = 1e-3
param = Parameter.defaultUnit(θ, rs)
kF, β = param.kF, param.β
Euv, rtol = 1000 * param.EF, 1e-11

# Nk, order = 8, 4
# maxK, minK = 10kF, 1e-7kF
# Nk, order = 11, 8
maxK, minK = 20kF, 1e-8kF
Nk, order = 16, 12
# maxK, minK = 10kF, 1e-9kF
# Euv, rtol = 1000 * param.EF, 1e-11
# maxK, minK = 20param.kF, 1e-9param.kF
# Nk, order = 16, 12

# test G0
kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [0.0, kF], Nk, minK, order)
G0 = SelfEnergy.G0wrapped(Euv, rtol, kgrid, param)
kF_label = searchsortedfirst(kgrid.grid, kF)
G_tau = GreenFunc.toTau(G0)
# println(G_tau.timeGrid.grid)
# println(real(G_tau.dynamic[1, 1, :, end]) .* (-1))
G_ins = tau2tau(G_tau.dlrGrid, G_tau.dynamic, [β,], G_tau.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
# println(real(G_ins))
integrand = real(G_ins) .* kgrid.grid .* kgrid.grid
density0 = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π^2
# println("$density0, $(param.n)")
# @test isapprox(param.n, density0, rtol=3e-5)

Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :ko)
Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

Z0 = (SelfEnergy.zfactor(Σ))
#z = zlist[ind]
# @test isapprox(Z0, z, rtol=3e-3)
mratio = SelfEnergy.massratio(param, Σ)
#m = mlist[ind]
# @test isapprox(mratio, m, rtol=3e-3)

println("θ = $θ,  rs= $rs")
println("Z-factor = $Z0")
println("m*/m = $mratio")
