
# benchmark with mc results

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

function oneloop(dim, θ, rs, channel; kwargs...)
    param = Parameter.defaultUnit(θ, rs, dim)
    if haskey(kwargs, :int_type) && kwargs[:int_type] == :none
        param = Parameter.Para(param; gs=0, ga=0)
    end
    uid = 0
    if haskey(kwargs, :uid)
        uid = kwargs[:uid]
    end

    println("dim=$dim, θ=$θ, rs=$rs, channel=$channel:")

    lamu, R_freq, F_freq, R_ins = BSeq.linearResponse(param, channel;
        Ntherm=0, Nmax=0, kwargs...)

    return R_freq, F_freq, R_ins
end

Rw, Fw, Ri = oneloop(3, 0.1, 3.0, 0; α=0.0)

println(Ri.mesh[2])
println(Ri.data)