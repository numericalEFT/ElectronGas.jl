
module MeasureChi

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

export measure_chi

function measure_chi(F_freq::GreenFunc.MeshArray)
    F_dlr = F_freq |> to_dlr
    F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
    integrand = view(F_ins, 1, :)
    kgrid = F_ins.mesh[2]
    return real(CompositeGrids.Interp.integrate1D(integrand, kgrid))
end

function measure_chi(dim, θ, rs)
    param = Parameter.rydbergUnit(θ, rs, dim)
    channel = 0
    lamu, R_freq, F_freq = BSeq.linearResponse(param, channel)
    return measure_chi(F_freq)
end

end

using Test
using .MeasureChi

@testset "measure chi" begin
    # println(measure_chi(3, 1e-2, 2.0))
    dim = 3
    rs = 3.0
    beta = [400, 800, 1600] .* 8
    chi = [measure_chi(dim, 1 / b, rs) for b in beta]
    println(chi)
end