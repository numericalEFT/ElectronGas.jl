using ElectronGas
using CompositeGrids

function Ws1(q, n, param)
    if abs(q) < 1e-16
        q = 1e-16
    end
    ωn = 2 * n * π / param.β
    return 2π * param.e0^2 / sqrt(q^2 + param.Λs) * ωn^2 / (ωn^2 + 4π * param.kF^2 * sqrt(q^2 + param.Λs))
end

function intWs2D(n, param; Ws=Ws1)
    xgrid = ElectronGas.CompositeGrids.CompositeG.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 24, 1e-12, 8)
    qs = [sqrt(param.kF^2 * 2 - 2 * cos(x) * param.kF^2) for x in xgrid]
    data = [Ws(q, n, param) for q in qs]
    return ElectronGas.CompositeGrids.Interp.integrate1D(data, xgrid)
end

param = ElectronGas.Parameter.defaultUnit(1e-5, 0.5, 2; Λs=1e-10)
# nlist = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
# nlist = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500]
nlist = [i for i in 0:10:5000]
wlist = [π * (2n + 1) / param.β for n in nlist]

W1 = ([intWs2D(n, param; Ws=Ws1) for n in nlist])
# println(W1)

using Plots
p = plot(wlist, W1)
display(p)
readline()
