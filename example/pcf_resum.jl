using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

using DelimitedFiles

function pcf_resum_ab(dim, θ, rs, channel; kwargs...)
    param = Parameter.rydbergUnit(θ, rs, dim; Λs=1e-6)
    if haskey(kwargs, :int_type) && kwargs[:int_type] == :none
        param = Parameter.Para(param; gs=0, ga=0)
    end
    if haskey(kwargs, :Ec_ratio)
        Ec = kwargs[:Ec_ratio] * param.EF
    end
    uid = 0
    if haskey(kwargs, :uid)
        uid = kwargs[:uid]
    end


    println("dim=$dim, θ=$θ, rs=$rs, channel=$channel:")

    # result = BSeq_resum.pcf_resum(param, channel; Ec=Ec, kwargs...)
    result = BSeq_resum.pcf_resum_smooth(param, channel; Ec=Ec, kwargs...)
    # dir = "./run/"
    # # fname = "gap$(dim)D_phchi_rs$(rs)_l$(channel)_vlargemu19.txt"
    # # fname = "gap$(dim)D_rpachi_rs$(rs)_l$(channel)_vcrit$(uid÷100).txt"
    # # fname = "gap$(dim)D_phrpachi_rs$(rs)_l$(channel)_vlarge0.txt"
    # # fname = "gap_plasmon_rs$(rs)_l$(channel)_vcrit$(uid÷100).txt"
    # fname = "gap_plasmonfs_rs$(rs)_l$(channel)_vcrit$(uid÷100).txt"
    # open(dir * fname, "a+") do io
    #     writedlm(io, data, ' ')
    # end

    # return 1 / result
    return result
end

using Test

@testset "pcf resum" begin
    # println(measure_chi(3, 1e-2, 2.0))
    uid0 = 3000087
    dim = 3
    rs = 3.0
    # num = 14
    # num = 25
    num = 1
    channel = 0
    # beta = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
    # beta = [400 * 2^(i - 1) for i in 1:num]
    beta = [6400,]
    # beta = [400 * 20000^(i / num) for i in LinRange(0, num - 1, num)]
    # beta = [400 * 20000^(i / num) for i in LinRange(0, num - 1, num)]
    # beta = [100 * 2^(i / num) for i in LinRange(0, num - 1, num)]
    # beta = [1.5625 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    # beta = [1.5625 * sqrt(2)^i for i in LinRange(17, num - 1, num - 17)]
    # beta = [800 * sqrt(2)^i for i in LinRange(0, 7 - 1, 7)]
    # beta = [1.5625 * 2^i for i in LinRange(0, num - 1, num)]
    # beta = [1800, 2000, 2229.78,]
    # num = 6
    # beta = [50 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    # chi = [measure_chi(dim, 1 / b, rs; sigmatype=:g0w0) for b in beta]
    result = [pcf_resum_ab(dim, 1 / beta[i], rs, channel;
        # atol=1e-8, rtol=1e-10, Nk=8, order=8, Ntherm=30, α=0.8,
        atol=1e-7, rtol=1e-10, Nk=8, order=4, Ntherm=5, α=0.75,
        # sigmatype=:none, int_type=:rpa, Vph=phonon,
        # sigmatype=:none, int_type=:rpa,
        sigmatype=:none, int_type=:ko,
        # sigmatype=:none, int_type=:ko, Vph=BSeq.phonon,
        # sigmatype=:none, int_type=:none, Vph=phonon,
        # plasmon_type=:plasmon,
        # plasmon_type=:plasmon_fs,
        # resum=true,
        Ec_ratio=6,
        onlyA=true,
        ω_c_ratio=0.1,
        issave=true, uid=uid0 + i, dir="./run/data/",
        verbose=true) for i in 1:length(beta)]
    # println(chi)
end
