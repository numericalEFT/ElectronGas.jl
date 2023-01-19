module Gamma

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids
using ElectronGas.Interaction
using DelimitedFiles

export measure_chi
# param = Parameter.rydbergUnit(100.0,1.0; γ = 0.5)
# γtest,γins = gamma_wrapped(100.0, 1e-10, param)

# γcompare = γtest |> to_dlr |>to_imfreq
# print("$(maximum(abs.(γtest.-γcompare)))")

# function measure_chi(F_freq::GreenFunc.MeshArray)
#     F_dlr = F_freq |> to_dlr
#     F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
#     kgrid = F_ins.mesh[2]
#     integrand = view(F_ins, 1, :)
#     return real(CompositeGrids.Interp.integrate1D(integrand, kgrid))
# end

function measure_chi( θ, g, γ; kwargs...)
    param = Parameter.defaultUnit(θ, 1.0; γ = γ, g = g)
    lamu, R_freq, F_freq = BSeq.linearResponse(param
        ; kwargs...)

    data = [1 / θ lamu]

    dir = "./run/"
    fname = "gap_gamma.txt"
    open(dir * fname, "a+") do io
        writedlm(io, data, ' ')
    end

    return lamu
end

end

using Test
using .Gamma

@testset "measure chi" begin
    # println(measure_chi(3, 1e-2, 2.0))
    dim = 3
    rs = 7.0
    num = 10
    #beta = [6.25 * 2^i for i in LinRange(0, num - 1, num)]
    # beta = [2000, 2200.0, 2229.78,]
    # beta = [200,]
    # num = 18
    beta = [200.0 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
    # chi = [measure_chi(dim, 1 / b, rs; sigmatype=:g0w0) for b in beta]
    g =-0.05/π
    γ = 1.0
    lamu = [measure_chi(1 / b, g, γ; atol=1e-6, rtol=1e-10, verbose=true, α=0.5) for b in beta]
    println(lamu)
end
