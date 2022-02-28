# Self-energy for spin-fermion model

using ElectronGas
using Test, Printf, DelimitedFiles
using Gaston, LaTeXStrings

bexp = 3
ExtQ = 0

dim = 2
beta, rs = 10^bexp, 1.0
ga = 1.0
param = Interaction.Parameter.rydbergUnit(1 / beta, rs, dim)

Λa = param.Λa + Polarization.Polarization0_ZeroTemp(0.0, 0, param) * param.spin * ga * (-param.e0a^2) / param.ϵ0
println(Λa)
param = Parameter.Para(param, gs = 0.0, ga = ga, Λa = Λa)

function sigma(param)
    Euv, rtol = 100 * param.EF, 1e-10
    Nk, order, minK = 8, 8, 1e-7
    # Nk, order, minK = 11, 4, 1e-8

    Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, minK * param.kF, order, :rpa)
    Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

    kgrid = Σ.spaceGrid
    println(kgrid)
    kF = kgrid.panel[3]
    kF_label = searchsortedfirst(kgrid.grid, kF)
    println(kF_label)
    ωgrid = Σ.dlrGrid

    ΣR = real(Σ.dynamic)
    ΣI = imag(Σ.dynamic)
    println(real(Σ.instant[1, 1, :]))
    println(ΣR[1, 1, kF_label, :])

    label = 1
    print([ωgrid.n ωgrid.ωn ΣR[1, 1, label, :] ΣI[1, 1, label, :]])
    # open("./plot/sigma_b1e$(bexp)_q$(ExtQ)kf_m0.txt", "w") do io
    #     writedlm(io, [ωgrid.n ωgrid.ωn ΣR[1, 1, label, :] ΣI[1, 1, label, :]])
    # end
end



if abspath(PROGRAM_FILE) == @__FILE__
    sigma(param)

    data = readdlm("./plot/sigma_b1e$(bexp)_q$(ExtQ)kf_m0.txt", '\t', Float64, '\n')
    # using Gaston
    # x = ωgrid.ωn / param.EF
    # y = ΣI[1, 1, kF_label, :]
    x = data[:, 2] / param.EF
    y = -data[:, end]
    plot(x, y, curveconf = "notit w lp lw 1 lc '#08F7FE'",
        Axes(xrange = (0, 0.1),
            # yrange = (-4, 1.2),
            ylabel = "'-Im Σ'",
            # ylabel = "'f_{xc} N_F (q/ω_n)^2'",
            xlabel = "'ω_n/E_F'",
            key = "t l",
            title = "'beta=$beta, r_s=$rs'")
    )
    x = 0:1e-2:2
    y = x .^ (2 / 3) .* 0.16
    plot!(x, y, curveconf = "tit 'ω_n^{2/3}'")
    save(term = "pdf", output = "./plot/SigmaIm_b1e$(bexp)_rs$(rs)_q$(ExtQ)kf.pdf", linewidth = 1)
end