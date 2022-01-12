@testset "Self Energy" begin
    param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(10000.0, 3.0)

    Euv, rtol = 100*param.EF, 1e-10
    kernel = SelfEnergy.LegendreInteraction.DCKernel(param, Euv, rtol, 10, 10*param.kF, 1e-7*param.kF, 5, :rpa,0,:sigma)

    G0 = SelfEnergy.G0wrapped(Euv, rtol, kernel.kgrid, param)

    Σ = SelfEnergy.calcΣ(G0, kernel)

    Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

    kgrid = Σ.spaceGrid
    kF = kgrid.panel[3]
    kF_label = searchsortedfirst(kgrid.grid, kF)

    ΣR = real(Σ.dynamic)
    ΣI = imag(Σ.dynamic)
    println(Σ.instant[1,1,:])
    println(ΣR[1,1,kF_label,:])
    println(ΣI[1,1,kF_label,:])
    println(SelfEnergy.zfactor(Σ))
end
