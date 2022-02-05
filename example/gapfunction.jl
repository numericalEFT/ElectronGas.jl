"""
Calculate gap-function equation
"""

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids

function G2wrapped(Σ::GreenFunc.Green2DLR, param)
    # return G(K)G(-K)
    @unpack me, kF, β, EF = param
    Σ_freq = GreenFunc.toMatFreq(Σ)
    Σ_shift = GreenFunc.dynamic(Σ_freq, π/β, kF, 1, 1) + GreenFunc.instant(Σ_freq, kF, 1, 1)
    green =  Green2DLR{ComplexF64}(
        :G, GreenFunc.IMFREQ,Σ_freq.β, Σ_freq.isFermi, Σ_freq.dlrGrid.Euv, Σ_freq.spaceGrid, Σ_freq.color;
        timeSymmetry = Σ_freq.timeSymmetry, rtol = Σ_freq.dlrGrid.rtol)

    green_dyn = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(green.spaceGrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            #green_dyn[1,1,ki,ni] = 1/(im*(π/β*(2n+1)) - (k^2/2/me-EF) + Σ.dynamic[1,1,ki,ni] + Σ.instant[1,1,ki])
            green_dyn[1,1,ki,ni] = 1/(
                ( (2n+1)*π/β-imag(Σ_freq.dynamic[1,1,ki,ni]) )^2
                + (ω + real(Σ_freq.dynamic[1,1,ki,ni] + Σ_freq.instant[1,1,ki] - Σ_shift) )^2
            )
        end
    end
    green.dynamic=green_dyn
    return green
end

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

function calcF(Δ::GreenFunc.Green2DLR, G2::GreenFunc.Green2DLR)
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
