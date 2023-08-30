using ElectronGas, JLD2
using ElectronGas.CompositeGrids
using ElectronGas.GreenFunc
using ElectronGas.Lehmann

fname = "run/data/sigma.jld2"
f = jldopen(fname, "r")
sigt, sigi = f["Dyn"], f["Ins"]
sigdlr = to_dlr(sigt)
sigw = to_imfreq(sigdlr)
wgrid = sigw.mesh[1]
kgrid = sigw.mesh[2]
ikF = searchsortedfirst(kgrid, 0.6397194308925043)
ΣI = imag(sigw.data[:, ikF] .+ sigi[1, ikF])
zph = [wgrid[i] / (wgrid[i] + ΣI[i]) for i in 1:length(wgrid)]
println(zph)

