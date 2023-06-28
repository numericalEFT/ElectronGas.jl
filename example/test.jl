using ElectronGas
using JLD2

function R0_freq(uid=3210005)
    fname = "run/data/PCFdata_$uid.jld2"
    f = jldopen(fname, "r")
    param = f["param"]
    ri, rw = f["R_ins"], f["R_freq"]
    kgrid = rw.mesh[2]
    kF = param.kF
    ikF = ElectronGas.GreenFunc.locate(kgrid, kF)
    # return ri[1, ikF] + real(rw[1, ikF]), ri[1, ikF] + real(rw[end, ikF])
    return ri[1, ikF] .+ real(rw[:, ikF])
end

fname = "run/data/PCFdata_3210005.jld2"
f = jldopen(fname, "r")
param = f["param"]
ri, rw = f["R_ins"], f["R_freq"]

wgrid = rw.mesh[1]
ngrid = wgrid.grid

