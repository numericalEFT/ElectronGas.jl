# Self-energy for spin-fermion model

using ElectronGas
using Test, Printf, DelimitedFiles

dim = 2
beta, rs = 1e2, 1.0
espin = 1.0
param = Interaction.Parameter.defaultUnit(1 / beta, rs, dim)

Λa = Polarization.Polarization0_ZeroTemp(0.0, 0, param) * (-espin^2) / param.ϵ0
println(Λa)

param = Parameter.Para(param, e0 = 0.0, espin = espin, Λa = Λa)

Euv, rtol = 100 * param.EF, 1e-10
# Nk, order, minK = 8, 4, 1e-7
Nk, order, minK = 11, 4, 1e-8

Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, minK * param.kF, order, :rpa)
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

# f = open("./data/Nk$Nk/SigmaIm_b1e2_m0.txt", "w")
for (n, sigma) in enumerate(ΣI[1, 1, kF_label, :])
    println(n, ' ', ωgrid.n[n], ' ', ωgrid.ωn[n], ' ', sigma)
    writedlm(f, [ωgrid.n[n] ωgrid.ωn[n] sigma])
end
# close(f)

println(SelfEnergy.zfactor(Σ))