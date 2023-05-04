# Self-energy for spin-fermion model

using ElectronGas, GreenFunc
using Test, Printf, DelimitedFiles
using Gaston, LaTeXStrings, Parameters

bexp = 3
ExtQ = 1

dim = 3
beta, rs = 10^bexp, 0.1
ga = 1.0
param = Interaction.Parameter.rydbergUnit(1 / beta, rs, dim)

Λa = param.Λa + Polarization.Polarization0_ZeroTemp(0.0, 0, param) * param.spin * ga * (-param.e0a^2) / param.ϵ0
println(Λa)
param = Parameter.Para(param, gs=0.0, ga=ga, Λa=Λa)

function sigma(param)
    @unpack me, kF, EF, β = param
    # Euv, rtol = 100 * param.EF, 1e-10
    Euv, rtol = 1000EF, 1e-11
    maxK, minK = 20kF, 1e-8kF
    Nk, order = 12, 8
    # Nk, order, minK = 8, 8, 1e-7
    # Nk, order, minK = 11, 4, 1e-8

    Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)
    # Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)
    Σ_freq = Σ |> to_dlr |> to_imfreq

    kgrid = Σ.mesh[2]
    # println(kgrid)
    kF_label = locate(kgrid, kF)
    k0_label = locate(kgrid, 0.0)
    println(kF_label)

    dlrgrid = Σ_freq.mesh[1].representation
    # println(dlrgrid.n)

    # println(Σ_freq[:, kF_label])
    # println(Σ_ins[1, kF_label])

    ΣR = real(Σ_freq[:, kF_label] .+ Σ_ins[1, kF_label])
    ΣI = imag(Σ_freq[:, kF_label])
    # println(real(Σ.instant[1, 1, :]))
    # println(ΣR[1, 1, kF_label, :])
    open("./plot/sigma$(dim)D_b1e$(bexp)_rs$(rs)_q$(ExtQ)kf_m0.txt", "w") do io
        writedlm(io, [dlrgrid.n dlrgrid.ωn ΣR ΣI])
        # writedlm(io, [ωgrid.n ωgrid.ωn ΣR[:, kF_label] ΣI[:, kF_label]])
    end

    kmax_label = locate(kgrid, 2.1 * kF)
    # kamp_grid = Vector{Float64}(undef, kmax_label)
    ds_dk = Vector{Float64}(undef, kmax_label)
    ds_dw = similar(ds_dk)
    Σ_freq = dlr_to_imfreq(to_dlr(Σ), [0, 1])
    sigma2 = real(Σ_freq[1, 1] + Σ_ins[1, 1])
    for (i, kamp) in enumerate(kgrid[1:kmax_label])
        k_label = locate(kgrid, kamp)
        sigma1 = real(Σ_freq[1, k_label] + Σ_ins[1, k_label])
        ds_dk[i] = (sigma1 - sigma2) / (kamp^2 / 2 / me)

        ΣI = imag(Σ_freq[:, k_label])
        ds_dw[i] = (ΣI[2] - ΣI[1]) / 2 / π * β
    end
    # label = 1
    open("./plot/Sigma$(dim)Dvsk_b1e$(bexp)_rs$(rs)_m0.txt", "w") do io
        writedlm(io, [kgrid[1:kmax_label] ds_dk ds_dw])
        # writedlm(io, [ωgrid.n ωgrid.ωn ΣR[:, kF_label] ΣI[:, kF_label]])
    end

    # open("./plot/sigma$(dim)D_b1e$(bexp)_q0kf_m0.txt", "w") do io
    #     writedlm(io, [ωgrid.n ωgrid.ωn ΣR[:, k0_label] ΣI[:, k0_label]])
    # end
end



if abspath(PROGRAM_FILE) == @__FILE__
    sigma(param)

    data = readdlm("./plot/sigma$(dim)D_b1e$(bexp)_rs$(rs)_q$(ExtQ)kf_m0.txt", '\t', Float64, '\n')
    # using Gaston
    x = data[:, 2] / param.EF
    y = -data[:, end]
    plot(x, y, curveconf="notit w lp lw 1 lc '#08F7FE'",
        Axes(xrange=(0, 1.0),
            # yrange = (-4, 1.2),
            ylabel="'-Im Σ(k_F,iω_n)'",
            # ylabel = "'f_{xc} N_F (q/ω_n)^2'",
            xlabel="'ω_n/E_F'",
            key="t l",
            title="'beta=$beta, r_s=$rs'")
    )
    # x = 0:1e-2:2
    # y = x .^ (2 / 3) .* 0.16
    # plot!(x, y, curveconf="tit 'ω_n^{2/3}'")
    save(term="pdf", output="./plot/Sigma$(dim)DIm_b1e$(bexp)_rs$(rs)_q$(ExtQ)kf.pdf", linewidth=1)

    data1 = readdlm("./plot/Sigma$(dim)Dvsk_b1e$(bexp)_rs$(rs)_m0.txt", '\t', Float64, '\n')
    x = data1[:, 1] / param.kF
    y1, y2 = data1[:, 2], data1[:, 3]
    plot(x, y1, curveconf="notit w lp lw 1 lc '#08F7FE'",
        Axes(yrange=(0, 0.3),
            ylabel="'(Σ_{k,iω_0}-Σ_{0,iω_0})/(k^2/2m)'",
            # ylabel = "'f_{xc} N_F (q/ω_n)^2'",
            xlabel="'k/k_F'",
            key="t l",
            title="'beta=$beta, r_s=$rs'")
    )
    save(term="pdf", output="./plot/dsdk$(dim)D_b1e$(bexp)_rs$(rs).pdf", linewidth=1)

    plot(x, y2, curveconf="notit w lp lw 1 lc '#08F7FE'",
        Axes(ylabel="'∂Σ_{k,iω_0}/∂ω'",
            # ylabel = "'f_{xc} N_F (q/ω_n)^2'",
            xlabel="'k/k_F'",
            key="t l",
            title="'beta=$beta, r_s=$rs'")
    )
    save(term="pdf", output="./plot/dsdw$(dim)D_b1e$(bexp)_rs$(rs).pdf", linewidth=1)
end