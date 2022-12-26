
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

function measure_chi(dim, θ, rs; kwargs...)
    param = Parameter.rydbergUnit(θ, rs, dim)
    channel = 0
    lamu, R_freq, F_freq = BSeq.linearResponse(param, channel
        ; kwargs...)
    result = measure_chi(F_freq)
    println("1/chi=", 1 / result)

    data = [1 / θ 1 / result lamu channel rs]

    dir = "./run/"
    fname = "gap_chi_rs$(rs)_l$(channel).txt"
    open(dir * fname, "a+") do io
        writedlm(io, data, ' ')
    end

    return 1 / result
end

end

using Test
using .MeasureChi

@testset "measure chi" begin
    # println(measure_chi(3, 1e-2, 2.0))
    dim = 3
    rs = 3.0
    num = 20
    beta = [6.25 * 2^i for i in LinRange(0, num - 1, num)]
    # beta = [2000, 2200.0, 2229.78,]
    # num = 18
    # beta = [6.25 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    # chi = [measure_chi(dim, 1 / b, rs; sigmatype=:g0w0) for b in beta]
    chi = [measure_chi(dim, 1 / b, rs; sigmatype=:none, atol=1e-6) for b in beta]
    println(chi)
end