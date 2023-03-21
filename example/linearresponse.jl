
module MeasureChi

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

using DelimitedFiles

export measure_chi

function measure_chi(F_freq::GreenFunc.MeshArray)
    F_dlr = F_freq |> to_dlr
    F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
    kgrid = F_ins.mesh[2]
    integrand = view(F_ins, 1, :)
    return real(CompositeGrids.Interp.integrate1D(integrand, kgrid))
end

function measure_chi(dim, θ, rs, channel; kwargs...)
    param = Parameter.rydbergUnit(θ, rs, dim)
    if haskey(kwargs, :int_type) && kwargs[:int_type] == :none
        param = Parameter.Para(param; gs=0, ga=0)
    end

    println("dim=$dim, θ=$θ, rs=$rs, channel=$channel:")

    lamu, R_freq, F_freq = BSeq.linearResponse(param, channel; kwargs...)
    result = measure_chi(F_freq)
    println("1/chi=", 1 / result)

    data = [1 / θ 1 / result lamu channel rs]

    dir = "./run/"
    fname = "gap$(dim)D_phchi_rs$(rs)_l$(channel)_vlargemu19.txt"
    # fname = "gap$(dim)D_rpachi_rs$(rs)_l$(channel)_v1.txt"
    # fname = "gap$(dim)D_phrpachi_rs$(rs)_l$(channel)_vlarge0.txt"
    open(dir * fname, "a+") do io
        writedlm(io, data, ' ')
    end

    return 1 / result
end

end

using Test
using .MeasureChi
using ElectronGas.Interaction

@testset "measure chi" begin
    # println(measure_chi(3, 1e-2, 2.0))
    dim = 3
    rs = 1.0
    # num = 14
    num = 25
    channel = 0
    # beta = [100 * 2^(i / num) for i in LinRange(0, num - 1, num)]
    # beta = [1.5625 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    beta = [1.5625 * sqrt(2)^i for i in LinRange(17, num - 1, num - 17)]
    # beta = [800 * sqrt(2)^i for i in LinRange(0, 7 - 1, 7)]
    # beta = [1.5625 * 2^i for i in LinRange(0, num - 1, num)]
    # beta = [1800, 2000, 2229.78,]
    # num = 6
    # beta = [50 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    # chi = [measure_chi(dim, 1 / b, rs; sigmatype=:g0w0) for b in beta]
    chi = [measure_chi(dim, 1 / b, rs, channel;
        atol=1e-8, rtol=1e-10, Nk=8, order=8, Ntherm=100, α=0.8,
        # sigmatype=:none, int_type=:rpa, Vph=phonon,
        # sigmatype=:none, int_type=:rpa,
        sigmatype=:none, int_type=:none, Vph=phonon,
        verbose=true) for b in beta]
    println(chi)
end