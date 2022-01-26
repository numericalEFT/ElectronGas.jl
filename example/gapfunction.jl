"""
Calculate gap-function equation
"""

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids



function calcΔ(F::GreenFunc.Green2DLR, W::LegendreInteraction.DCKernel)
    @unpack β = W.param
    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = F.dlrGrid
    bdlr = W.dlrGrid

    Δ = GreenFunc.Green2DLR{Float64}(:delta, GreenFunc.IMTIME,β,true,fdlr.Euv,kgrid,1)
    Δ.instant = ones(Float64, (1,1,length(kgrid.grid)))
    Δ.dynamic = ones(Float64, (1,1,length(kgrid.grid), fdlr.size))

    return Δ
end

function calcF(Δ::GreenFunc.Green2DLR, G::GreenFunc.Green2DLR)
    return Δ
end

function gapIteration(param, N, Euv, rtol, Nk, maxK, minK, order, int_type)
    @unpack β = param

    W = SelfEnergy.LegendreInteraction.DCKernel0(param;
                                                      Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma)
    G = SelfEnergy.G0wrapped(Euv, rtol, W.kgrid, param)

    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = G.dlrGrid
    bdlr = W.dlrGrid

    Δ = GreenFunc.Green2DLR{Float64}(:delta, GreenFunc.IMTIME,param.β,true,G.dlrGrid.Euv,W.kgrid,1)
    Δ.instant = ones(Float64, (1,1,length(kgrid.grid)))
    Δ.dynamic = ones(Float64, (1,1,length(kgrid.grid), fdlr.size))

    return Δ
end

if abspath(PROGRAM_FILE) == @__FILE__
    param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(1/1000.0, 1.0)
    Euv, rtol = 100*param.EF, 1e-10
    Nk, order = 8, 4

    Δ = gapIteration(param, 1, Euv, rtol, Nk, 10*param.kF, 1e-7*param.kF, order, :rpa)

    kgrid = Δ.spaceGrid
    kF = kgrid.panel[3]
    kF_label = searchsortedfirst(kgrid.grid, kF)

    println(Δ.instant[1,1,:])
    println(Δ.dynamic[1,1,kF_label,:])
end
