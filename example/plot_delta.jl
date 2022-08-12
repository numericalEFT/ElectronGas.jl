using LinearAlgebra
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids
using ElectronGas.Convention
using DelimitedFiles
using Plots
using JLD2, FileIO
using CodecZlib
using LaTeXStrings

dir = "./run/"
# fnames = readdir(dir)
# for fname in fnames
#     if fname[end-3:end] == "jld2"
#         println(fname)
#     end
# end

fname = "gap_minisub_rs2.2_l0.jld2"
d = load(dir * fname)


deltas = d["deltas"]

plt = plot(xlims = (0.0, 20.0), ylims = (-7, 1))
# plt = plot()
# xlabel!(plt, L"$ω_n/E_F\ln{(\beta T_F)}$")
# ylabel!(plt, L"$[\Delta(ω_n, k_F) - \Delta(ω1,k_F)]/ \ln{(\beta T_F)}^2$")
xlabel!(plt, L"$ω_n/E_F$")
ylabel!(plt, L"$\Delta(ω_n, k_F)$")

betas = zeros(length(deltas))
deltainfs = zeros(length(deltas))

for (i, delta) in enumerate(deltas)
    kgrid = delta.spaceGrid
    kF = 1.0
    kF_label = searchsortedfirst(kgrid.grid, kF)
    omegagrid = delta.dlrGrid.ωn
    beta = delta.dlrGrid.β
    deltakF = real(delta.dynamic[1, 1, kF_label, :] .+ delta.instant[1, 1, kF_label])
    # nmax = 15
    # plot!(plt, omegagrid[1:nmax] .* log(beta), deltakF[1:nmax], label = L"$\beta T_F=%$beta$")
    # plot!(plt, omegagrid .* (log(beta)), (deltakF .- deltakF[1]) ./ (log(beta))^2, label = L"$\beta T_F=%$beta$")
    plot!(plt, omegagrid, deltakF, label=L"$\beta T_F=%$beta$")
    betas[i], deltainfs[i] = beta, deltakF[end]
end

println(betas)
println(deltainfs)

display(plt)
savefig(plt, "delta_log.pdf")

plt2 = plot()
plot!(plt2, log.(betas) .^ 2, deltainfs)
display(plt2)
readline()
