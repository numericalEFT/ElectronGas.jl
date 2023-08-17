using ElectronGas, JLD2
using CompositeGrids
using LsqFit

function D_freq(ωn, param)
    ωD = wdratio * param.EF
    # ωn = 2 * n * π / param.β
    return ωD^2 / (ωn^2 + ωD^2)
end

f = jldopen("run/data/PCFdata_3500003.jld2", "r")
# f = jldopen("run/data/PCFdata_300005.jld2", "r")

param = f["param"]
ri = f["R_ins"]
rw = f["R_freq"]
kF = param.kF

println("β=$(param.β)")

kgrid = rw.mesh[2]

iw0 = 10
rw0 = ri.data[1, :] .+ rw.data[iw0, :]
println(rw0)
ikF = searchsortedfirst(kgrid, param.kF)
println(ikF)

@. model(x, p) = p[1] + p[2] * (x - kF)^2 + p[3] * (x - kF)^4
dk = 0.2
kfgrid = SimpleG.Uniform([param.kF * (1 - dk), param.kF * (1 + dk)], 100)
data = Interp.interp1DGrid(rw0, kgrid, kfgrid)

fit = curve_fit(model, kfgrid, data, [1.0, 0.0, 0.0])
println(coef(fit))

using Plots
# p = plot(kgrid[ikF-dn:ikF+dn], rw0[ikF-dn:ikF+dn])
p = plot(kfgrid, data)
plot!([kF, kF], [0.0, 0.1])
display(p)
readline()