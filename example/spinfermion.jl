# Self-energy for spin-fermion model

using ElectronGas
using Test, Printf, DelimitedFiles

β = 1e3
rs = 1.0
espin = 1.0
# param = Interaction.Parameter.defaultUnit(beta, rs)
param = Interaction.Parameter.rydbergUnit(β, rs)
# param_b = Parameter.Para(param, e0 = 0.0, espin = 1.0)

Λa = Polarization.Polarization0_ZeroTemp(0.0, 0, param) * (espin^2)
println(Λa)

param = Parameter.Para(param, e0 = 0.0, espin = espin, Λa = Λa)

Euv, rtol = 100 * param.EF, 1e-10
Nk, order = 8, 4

Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, 1e-7 * param.kF, order, :rpa)
Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

kgrid = Σ.spaceGrid
kF = kgrid.panel[3]
kF_label = searchsortedfirst(kgrid.grid, kF)
println(kF_label)
ωgrid = Σ.dlrGrid

ΣR = real(Σ.dynamic)
ΣI = imag(Σ.dynamic)
println(Σ.instant[1, 1, :])
println(ΣR[1, 1, kF_label, :])

f = open("./data/SigmaIm.txt", "w")
for (n, ω) in enumerate(ΣI[1, 1, kF_label, :])
    # if ωgrid.n[n] < 0
    #     continue
    # end
    println(n, ' ', ωgrid.n[n], ' ', ωgrid.ωn[n], ' ', ω)
    writedlm(f, [ωgrid.n[n] ωgrid.ωn[n] ω])
end
close(f)

println(SelfEnergy.zfactor(Σ))