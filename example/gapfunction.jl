
using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

using DelimitedFiles

function gapfunction(dim, θ, rs, channel; kwargs...)
    param = Parameter.rydbergUnit(θ, rs, dim)
    if haskey(kwargs, :int_type) && kwargs[:int_type] == :none
        param = Parameter.Para(param; gs=0, ga=0)
    end
    uid = 0
    if haskey(kwargs, :uid)
        uid = kwargs[:uid]
    end


    println("dim=$dim, θ=$θ, rs=$rs, channel=$channel:")

    lamu, R_freq, F_freq = BSeq.lin_eliashberg(param, channel; kwargs...)

    data = [1 / θ lamu channel rs]

    dir = "./run/"
    # fname = "gap$(dim)D_phchi_rs$(rs)_l$(channel)_vlargemu19.txt"
    fname = "eigen$(dim)D_rpachi_rs$(rs)_l$(channel)_vcrit$(uid÷100).txt"
    # fname = "gap$(dim)D_phrpachi_rs$(rs)_l$(channel)_vlarge0.txt"
    open(dir * fname, "a+") do io
        writedlm(io, data, ' ')
    end

    return lamu
end

uid0 = 18500
dim = 3
rs = 1.85
# num = 14
# num = 25
num = 12
channel = 0
# beta = [400, 800]
beta = [400 * 2^(i - 1) for i in 1:num]
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
lamus = [gapfunction(dim, 1 / beta[i], rs, channel;
    shift=3.0,
    atol=1e-9, rtol=1e-10, Nk=8, order=8, Ntherm=1, α=0.8,
    # sigmatype=:none, int_type=:rpa, Vph=phonon,
    sigmatype=:none, int_type=:rpa,
    # sigmatype=:none, int_type=:none, Vph=phonon,
    issave=true, uid=uid0 + i, dir="./run/data/",
    verbose=true) for i in 1:length(beta)]
println(lamus)
