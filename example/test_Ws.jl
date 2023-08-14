using ElectronGas, JLD2
using ElectronGas.CompositeGrids

function Ws1(q, n, param)
    if abs(q) < 1e-16
        q = 1e-16
    end
    ωn = 2 * n * π / param.β
    return 4π * param.e0^2 / (q^2 + param.Λs) * ωn^2 / (ωn^2 + param.ωp^2 / 2)
end

@inline function Vinstant(q, param)
    if abs(q) < 1e-16
        q = 1e-16
    end
    return 4π * param.e0^2 / (q^2 + param.Λs)
end

@inline function Vcoulomb(q, param)
    if abs(q) < 1e-16
        q = 1e-16
    end
    return 4π * param.e0^2 / q^2
end

function KO_Ws_reg(n::Integer, para)
    q = 1e-16
    Pi = Polarization.Polarization0_3dZeroTemp(q, n, para)
    V = Vcoulomb(q, para)
    return 1 / (1 - Pi * V)
end

function Ws2(q, n, para)
    return KO_Ws_reg(n, para) * Vinstant(q, para)
end

function FSavg(n, para; Ws=Ws2)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 32, 1e-12, 8)
    qs = [sqrt(kF^2 * 2 - 2 * x * kF^2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = Ws(q, n, para)
    end
    return Interp.integrate1D(Wp, xgrid)
end

function intWs(n, param; Ws=Ws1)
    xgrid = ElectronGas.CompositeGrids.CompositeG.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 32, 1e-16, 8)
    qs = [sqrt(param.kF^2 * 2 - 2 * x * param.kF^2) for x in xgrid]
    data = [Ws(q, n, param) for q in qs]
    return ElectronGas.CompositeGrids.Interp.integrate1D(data, xgrid)
end

# f = jldopen("run/data/PCFresumdlr_300600001.jld2", "r")
# f = jldopen("run/data/PCFresumdlr_30000200001.jld2", "r")
f = jldopen("run/data/PCFresumdlr_3000200001.jld2", "r")
param, B = f["param"], f["B"]
nlist = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
wlist = [π * (2n + 1) / param.β for n in nlist]

data = B[1, :, 1]
wgrid = B.mesh[1]
println(data)

ω_c = 0.1param.EF
factor = (1 / 4 / π^2 * log(ω_c * param.β / 0.882))
ΠB00 = B[1] * param.EF * factor
println(ΠB00)

Wb = ([Interp.interp1D(data, wgrid, w) for w in wlist])
W1 = ([intWs(n, param; Ws=Ws1) for n in nlist])
W2 = ([intWs(n, param; Ws=Ws2) for n in nlist])
println(Wb)
println(W1)
println(W2)
# println(Wb .+ W2)
println(Wb .+ W2 .* (1 - ΠB00))