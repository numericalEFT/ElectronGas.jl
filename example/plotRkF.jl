using Plots
using JLD2
using ElectronGas
using ElectronGas.CompositeGrids
using ElectronGas.GreenFunc
using Roots

function load_PCFdata(fname)
    data = jldopen(fname, "r")

    param = data["param"]
    lamu = data["lamu"]
    F_freq = data["F_freq"]
    R_ins = data["R_ins"]
    R_freq = data["R_freq"]

    return param, lamu, F_freq, R_ins, R_freq
end

function extract_Rkf(param, R_ins, R_freq)
    ωgrid = R_freq.mesh[1]
    kgrid = R_freq.mesh[2]

    kF = param.kF
    ikF = searchsortedfirst(kgrid, kF)

    Rkf_ins = R_ins[ikF]
    Rkf_freq = R_freq[:, ikF]

    Rkf = real(Rkf_freq) .+ Rkf_ins
    return ωgrid, Rkf
end

function Nfind_zero(data, grid)
    wgrid = CompositeGrids.SimpleG.Arbitrary([w for w in grid])
    f(x) = CompositeGrids.Interp.interp1D(data, wgrid, x)
    w0 = find_zero(f, (wgrid[1], wgrid[end]))
    return w0
end

function measure_chi(F_freq::GreenFunc.MeshArray)
    F_dlr = F_freq |> to_dlr
    F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
    kgrid = F_ins.mesh[2]
    integrand = view(F_ins, 1, :)
    return real(CompositeGrids.Interp.integrate1D(integrand, kgrid))
end

uid0 = 2000000
Nrun = 12
params, lamus, Fs, Ris, Rfs = [], [], [], [], []
Rkfs = []
for i in 1:Nrun
    uid = uid0 + i
    fname = "run/data/PCFdata_$(uid).jld2"
    param, lamu, F_freq, R_ins, R_freq = load_PCFdata(fname)
    push!(params, param)
    push!(lamus, lamu)
    push!(Fs, F_freq)
    push!(Ris, R_ins)
    push!(Rfs, R_freq)
    wg, Rk = extract_Rkf(param, R_ins, R_freq)
    push!(Rkfs, Rk)
end

lnbetas = [log10(param.beta) for param in params]
Rinfs = [R[end] for R in Rkfs]
w0s = [Nfind_zero(Rkfs[i], Rfs[i].mesh[1]) for i in 1:Nrun]
chis = [measure_chi(Fs[i]) for i in 1:Nrun]
println(lnbetas)
println(Rinfs)
println(w0s)
println(chis)

p = plot(lnbetas, w0s .* Rinfs)
# p = plot(lnbetas, (chis))
# plot!(lnbetas, lamus)
display(p)
readline()

# fname = "data/PCFdata_1808.jld2"

# param, lamu, F_freq, R_ins, R_freq = load_PCFdata(fname)

# ωgrid = R_freq.mesh[1]
# kgrid = R_freq.mesh[2]

# kF = param.kF
# ikF = searchsortedfirst(kgrid, kF)

# Rkf_ins = R_ins[ikF]
# Rkf_freq = R_freq[:, ikF]

# Rkf = real(Rkf_freq) .+ Rkf_ins
# Fkf = real(F_freq[:, ikF])

# Rhalf = (Rkf[1] + Rkf[end]) / 2
# wgrid = CompositeGrids.SimpleG.Arbitrary([w for w in ωgrid])
# f(x) = CompositeGrids.Interp.interp1D(Rkf, wgrid, x)
# w0 = find_zero(f, (wgrid[1], wgrid[end]))
# wh = find_zero(x -> f(x) - Rhalf, (wgrid[1], wgrid[end]))
# println("rs=$(param.rs), β=$(param.beta)")
# println("R0, R∞ = ($(Rkf[1]), $(Rkf[end])), R0/R∞ = $(Rkf[1]/Rkf[end])")
# println("w0=$w0")
# println("wh=$wh")

# p = plot(ωgrid, Rkf)
# display(p)
# readline()