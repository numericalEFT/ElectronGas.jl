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
function AB_vals(uid)
    fname = "run/data/PCFresumdlr_$uid.jld2"
    f = jldopen(fname, "r")
    param = f["param"]
    A = f["A"]
    B = f["B"]
    wgrid = A.mesh[1]
    return param, wgrid, A.data, B
end

using ElectronGas.GreenFunc
using ElectronGas.GreenFunc.CompositeGrids
Πs(ω; ω_c=0.1) = atan(ω_c / abs(ω)) / abs(ω)
β = 3200
α = 0.882
Nmax = 100000
wg1 = GreenFunc.ImFreq(β, FERMION; grid=[i for i in 0:Nmax])
wg2 = CompositeG.LogDensedGrid(:cheb, [α / β, π / β * (2 * Nmax + 1)], [α / β, 0.1], 5, 0.001, 5)
data1 = [Πs(wg1[i]) for i in 1:length(wg1)]
data2 = [Πs(wg2[i]) for i in 1:length(wg2)]
int1 = sum(data1) / β
int2 = Interp.integrate1D(data2, wg2) / 2π
println((int1, int2, log(0.1 * β) / π))
# using Plots
# p = plot(xlim=(0.0, 5.0))
# for duid in 1:5
#     uid = duid + 3000000
#     param, wg, A, B = AB_vals(uid)
#     println(B)
#     plot!(p, wg, A)
# end
# display(p)
# readline()

# for duid in 1:9
#     uid = duid + 1230000
#     A, w0 = A_vals(uid)
#     println((A[1], A[end], w0))
# end

# using Plots
# p = plot(xlim=(0.0, 10.0))
# for duid in 1:9
#     uid = duid + 1230000
#     A, w0 = A_vals(uid)
#     println((A[1], A[end], w0))
#     plot!(p, A.mesh[1], A.data)
# end
# display(p)
# readline()


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

# using ElectronGas, JLD2
# f = jldopen("run/data/PCFresumdlr_3000010.jld2", "r")
# f["B"].mesh[1].grid[1:3:end]'
# f["B"].data[1:3:end, 1]'
