using ElectronGas
using JLD2
function A_vals(uid=3210005)
    fname = "run/data/PCFresum_$uid.jld2"
    f = jldopen(fname, "r")
    param = f["param"]
    A = f["A"]
    wgrid = A.mesh[1]
    iw0 = findfirst([(A[i]) <= 0.0 for i in 1:length(A)])
    return A, wgrid[iw0] / param.EF
end
# for duid in 1:9
#     uid = duid + 1230000
#     A, w0 = A_vals(uid)
#     println((A[1], A[end], w0))
# end

using Plots
p = plot(xlim=(0.0, 10.0))
for duid in 1:9
    uid = duid + 1230000
    A, w0 = A_vals(uid)
    println((A[1], A[end], w0))
    plot!(p, A.mesh[1], A.data)
end
display(p)
readline()


# function R0_freq(uid=3210005)
#     fname = "run/data/PCFdata_$uid.jld2"
#     f = jldopen(fname, "r")
#     param = f["param"]
#     ri, rw = f["R_ins"], f["R_freq"]
#     kgrid = rw.mesh[2]
#     kF = param.kF
#     ikF = ElectronGas.GreenFunc.locate(kgrid, kF)
#     # return ri[1, ikF] + real(rw[1, ikF]), ri[1, ikF] + real(rw[end, ikF])
#     return ri[1, ikF] .+ real(rw[:, ikF])
# end

# fname = "run/data/PCFdata_3210005.jld2"
# f = jldopen(fname, "r")
# param = f["param"]
# ri, rw = f["R_ins"], f["R_freq"]

# wgrid = rw.mesh[1]
# ngrid = wgrid.grid
