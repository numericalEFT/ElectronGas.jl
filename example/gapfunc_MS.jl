"""
Calculate gap-function equation
"""

using LinearAlgebra
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids
using ElectronGas.Convention
using DelimitedFiles

include("./eigen.jl")

const dim = 3
const freq_sep = 0.01

function G02wrapped(Euv, rtol, sgrid, param)
    # return G0(K)G0(-K)
    @unpack me, kF, β, μ = param

    green = GreenFunc.Green2DLR{Float64}(:g02, GreenFunc.IMFREQ, β, true, Euv, sgrid, 1; timeSymmetry=:pha, rtol=rtol)
    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, ωn) in enumerate(green.dlrGrid.ωn)
            ω = k^2 / 2 / me - μ
            green_dyn[1, 1, ki, ni] = 1 / (ωn^2 + ω^2)
        end
    end
    green.dynamic = green_dyn
    return green
end

function G2wrapped(Σ::GreenFunc.Green2DLR, param)
    # return G(K)G(-K)
    @unpack me, kF, β, μ = param

    green = Green2DLR{Float64}(
        :g2, GreenFunc.IMFREQ, β, true, Σ.dlrGrid.Euv, Σ.spaceGrid, 1;
        timeSymmetry=:pha, rtol=Σ.dlrGrid.rtol)
    Σ_freq = GreenFunc.toMatFreq(Σ, green.dlrGrid.n)

    green_dyn = zeros(Float64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    Σ_shift = real(GreenFunc.dynamic(Σ_freq, π / β, kF, 1, 1) + GreenFunc.instant(Σ_freq, kF, 1, 1))

    for (ki, k) in enumerate(green.spaceGrid)
        for (ni, ωn) in enumerate(green.dlrGrid.ωn)
            ω = k^2 / 2 / me - μ
            ΣR, ΣI = real(Σ_freq.dynamic[1, 1, ki, ni] + Σ_freq.instant[1, 1, ki] - Σ_shift), imag(Σ_freq.dynamic[1, 1, ki, ni])
            green_dyn[1, 1, ki, ni] = 1 / ((ωn - ΣI)^2 + (ω + ΣR)^2)
        end
    end
    green.dynamic = green_dyn
    return green
end

function ΔFinit(Euv, rtol, sgrid, param)
    @unpack β, me, μ = param
    Ω_c = freq_sep

    Δ = GreenFunc.Green2DLR{Float64}(:delta, GreenFunc.IMTIME, β, true, Euv, sgrid, 1; timeSymmetry=:pha, rtol=rtol)
    F = GreenFunc.Green2DLR{Float64}(:F, GreenFunc.IMFREQ, β, true, Euv, sgrid, 1; timeSymmetry=:pha, rtol=rtol)

    Δ_ins = ones(Float64, (Δ.color, Δ.color, Δ.spaceGrid.size))
    Δ_dyn = ones(Float64, (Δ.color, Δ.color, Δ.spaceGrid.size, Δ.timeGrid.size))
    F_dyn = zeros(Float64, (F.color, F.color, F.spaceGrid.size, F.timeGrid.size))

    for (ki, k) in enumerate(Δ.spaceGrid.grid)
        for (ni, ωn) in enumerate(Δ.dlrGrid.ωn)
            e = k^2 / 2 / me - μ
            Δ_dyn[1, 1, ki, ni] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
        end
    end
    Δ_dyn = real(matfreq2tau(Δ.dlrGrid, Δ_dyn, Δ.dlrGrid.τ, F.timeGrid.grid; axis=4))

    Δ.dynamic, Δ.instant = Δ_dyn, Δ_ins
    F.dynamic = F_dyn
    return Δ, F
end

function dotΔ(Δ, Δ0, Δ2=Δ, Δ02=Δ0)
    # kF_label = searchsortedfirst(kgrid.grid, param.kF)
    # return Δ0[kF_label] + real(Lehmann.tau2matfreq(fdlr, view(Δ, kF_label, :), fdlr.n))[1]
    # return (dot(Δ, Δ2) + dot(Δ0, Δ02))
    return (dot(Δ, Δ2))
end

function gapIteration_Renorm(param, channel, G2::GreenFunc.Green2DLR, kernel, kernel_bare, qgrids, Euv;
    Ntherm=200, rtol=1e-10)
    @unpack dim, kF = param

    kgrid = G2.spaceGrid

    Δ, F = ΔFinit(Euv, G2.dlrGrid.rtol, kgrid, param)
    delta = zeros(Float64, (kgrid.size, Δ.dlrGrid.size))
    delta0 = zeros(Float64, (kgrid.size))

    Δ0, Δ0_sum = 1.0, 0.0
    dΔ0, dΔ0_sum = zeros(Float64, (kgrid.size)), zeros(Float64, (kgrid.size))
    Δdyn_sum = zeros(Float64, (kgrid.size, Δ.dlrGrid.size))

    n = 0
    lamu, lamu0 = 1.0, 2.0
    err = 1.0
    kF_label = searchsortedfirst(kgrid.grid, kF)
    α = 0.8

    printn = 10
    while (true)
        n = n + 1
        delta = copy(Δ.dynamic[1, 1, :, :])
        delta0 = copy(Δ.instant[1, 1, :])

        calcF!(F, Δ, G2)
        if dim == 3
            calcΔ!(F, Δ, kernel, kernel_bare, qgrids)
        elseif dim == 2
            calcΔ_2d!(F, Δ, kernel, kernel_bare, qgrids)
        end

        Δ_kF = real(tau2matfreq(Δ.dlrGrid, Δ.dynamic, [0,], Δ.timeGrid.grid; axis=4)[1, 1, kF_label, 1] + Δ.instant[1, 1, kF_label])

        # if n < Ntherm
        Δ0_sum = kF + Δ_kF + Δ0_sum * α
        dΔ0_sum = Δ.instant[1, 1, :] + kgrid.grid .- (kF + Δ_kF) + dΔ0_sum .* α
        Δ0 = Δ0_sum * (1 - α)
        dΔ0 = dΔ0_sum .* (1 - α)
        Δ.instant[1, 1, :] = dΔ0 .+ Δ0

        Δdyn_sum = Δ.dynamic[1, 1, :, :] + Δdyn_sum .* α
        Δ.dynamic[1, 1, :, :] = Δdyn_sum .* (1 - α)
        # else
        #     Δ0 = kF + Δ_kF
        #     Δdyn_sum = Δ.dynamic[1, 1, :, :] + Δdyn_sum .* α
        #     Δ.dynamic[1, 1, :, :] = Δdyn_sum .* (1 - α)
        # end
        println("Δ(kF, ω0) = $Δ_kF, Δ0 = $Δ0  ($(kF + Δ_kF))")
        # n > Ntherm && lamu > lamu0 > 0 && break
        # n > Ntherm && lamu < lamu0 < 0 && break
        if n > Ntherm
            lamu = 1 - 1 / (1 + Δ_kF / kF)
            err = abs(lamu - lamu0) / abs(lamu + EPS)
            # lamu > lamu0 > 0 && err < rtol && break
            err < rtol && break
            lamu0 = lamu
        end
        if n%printn == 0
            println(1 - kF / Δ0)
        end
    end
    println("α = $α, iteration step: $n")
    lamu = 1 - kF / Δ0
    return lamu, Δ, F
end

function gapIteration(param, G2, kernel, kernel_bare, qgrids, Euv;
    Nstep=1e3, rtol=1e-6, shift=2.0)
    @unpack dim = param

    Δ, F = ΔFinit(Euv, G2.dlrGrid.rtol, G2.spaceGrid, param)
    delta = zeros(Float64, (G2.spaceGrid.size, Δ.dlrGrid.size))
    delta_sum = zeros(Float64, (G2.spaceGrid.size, Δ.dlrGrid.size))
    delta0 = zeros(Float64, (G2.spaceGrid.size))
    delta0_sum = zeros(Float64, (G2.spaceGrid.size))

    n = 0
    lamu, lamu0 = 1.0, 2.0
    err = 1.0

    while (n < Nstep && err > rtol)
        n = n + 1
        delta = copy(Δ.dynamic[1, 1, :, :])
        delta0 = copy(Δ.instant[1, 1, :])

        calcF!(F, Δ, G2)
        if dim == 3
            calcΔ!(F, Δ, kernel, kernel_bare, qgrids)
        elseif dim == 2
            calcΔ_2d!(F, Δ, kernel, kernel_bare, qgrids)
        end

        lamu = dotΔ(Δ.dynamic[1, 1, :, :], Δ.instant[1, 1, :], delta, delta0)

        Δ.dynamic[1, 1, :, :] = Δ.dynamic[1, 1, :, :] .+ shift .* delta
        Δ.instant[1, 1, :] = Δ.instant[1, 1, :] .+ shift .* delta0

        modulus = sqrt(dotΔ(Δ.dynamic[1, 1, :, :], Δ.instant[1, 1, :]))
        Δ.dynamic = Δ.dynamic ./ modulus
        Δ.instant = Δ.instant ./ modulus

        err = abs(lamu - lamu0) / abs(lamu + EPS)
        lamu0 = lamu
        # println(lamu)
    end
    return lamu, Δ, F
end


function gapfunction(beta, rs, channel::Int, dim::Int; sigmatype=:none, methodtype=:explicit, Ntherm = 20, int_type=:rpa, Λs=Λs)
    #--- parameters ---
    param = Parameter.defaultUnit(1 / beta, rs, dim; Λs=Λs)
    # param = Parameter.rydbergUnit(1 / beta, rs, dim)
    # Euv, rtol = 100 * param.EF, 1e-11
    # maxK, minK = 20param.kF, 1e-8param.kF
    # Nk, order = 12, 10
    Euv, rtol = 100 * param.EF, 1e-10
    maxK, minK = 10param.kF, 1e-7param.kF
    Nk, order = 10, 6

    #--- prepare kernel ---
    if dim == 3
        if channel == 0
            @time W = LegendreInteraction.DCKernel0(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
                int_type=int_type, channel=channel)
        else
            @time W = LegendreInteraction.DCKernel_old(param;
                Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
                int_type=int_type, channel=channel)
        end
    elseif dim == 2
        @time W = LegendreInteraction.DCKernel_2d(param;
            Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
            int_type=int_type, channel=channel)
    else
        error("No support for $dim dimension!")
    end

    fdlr = Lehmann.DLRGrid(Euv, param.β, rtol, true, :pha)
    bdlr = W.dlrGrid
    kgrid = W.kgrid
    qgrids = W.qgrids

    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3))

    kF_label = searchsortedfirst(kgrid.grid, param.kF)
    qF_label = searchsortedfirst(qgrids[kF_label].grid, param.kF)

    println("static kernel at (kF, kF):$(kernel_bare[kF_label, qF_label])")
    # println("dynamic kernel at (kF, kF):")
    # println(view(kernel_freq, kF_label, qF_label, :))

    #--- prepare Σ and G2 ---
    if sigmatype == :none
        G2 = G02wrapped(Euv, rtol, kgrid, param)
    elseif sigmatype == :g0w0
        Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        G2 = G2wrapped(Σ, param)
    end
    println("G2 at kF:")
    println(view(G2.dynamic,1, 1, kF_label, :))

    if methodtype == :explicit
        lamu, Δ, F = gapIteration(param, G2, kernel, kernel_bare, qgrids, Euv; rtol=1e-7, shift=3.0)
    elseif methodtype == :minisub
        lamu, Δ, F = gapIteration_Renorm(param, channel, G2, kernel, kernel_bare, qgrids, Euv;
                                         rtol=rtol, Ntherm=Ntherm)
    else
        error("method $(methodtype) not implemented!")
    end

    println("lamu = $lamu")

    Δ_freq = GreenFunc.toMatFreq(Δ)
    println(fdlr.ωn ./ param.EF)
    println(view(real(Δ_freq.dynamic), 1, 1, kF_label, :))
    # println(view(F, kF_label, :))

    data = [beta lamu channel rs]

    dir = "./run/"
    fname = "gap_$(methodtype)_rs$(rs)_l$(channel).txt"
    open(dir * fname, "a+") do io
        writedlm(io, data, ' ')
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # sigmatype = :none
    methodtype = :explicit
    sigmatype = :g0w0
    int_type = :rpa

    rs = 1.5
    channel = 0
    # channels = [0, 0, 0]
    # channels = [0, 1, 2, 3]

    if !isempty(ARGS)
        rs = parse(Float64, ARGS[1])
        channel = parse(Int, ARGS[2]) - 1
        if length(ARGS) == 3
            Λs = parse(Float64, ARGS[3])
        end
    end

    num = 9
    # blist = [400, 800, 1600, 3200, 6400]
    # blist = [6.25 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    blist = [400 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    # blist = [400 * 2^i for i in LinRange(0, num - 1, num)]
    # blist = [1 / 1.89059095e-05, 1 / 8.44687571e-05, 1 / 1.28551713e-05, 1 / 1.14498145e-06]
    # blist = [1 / 1.78760981e-05, 1 / 8.35387289e-05, 1 / 1.20811357e-05, 1 / 1.03562748e-06]

    for (ind, beta) in enumerate(blist)
    # for (channel, beta) in zip(channels, blist)
        println("beta = $beta,    rs = $rs")
        println("channel = $channel")
        gapfunction(beta, rs, channel, dim; sigmatype=sigmatype, methodtype=methodtype, int_type = int_type, Λs = Λs)
    end
end
