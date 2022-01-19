@testset "Self Energy" begin
    param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(1000.0, 1.0)
    Euv, rtol = 100 * param.EF, 1e-10
    Nk, order = 8, 4

    # kernel = SelfEnergy.LegendreInteraction.DCKernel(param, Euv, rtol, 10, 10*param.kF, 1e-7*param.kF, 5, :rpa,0,:sigma)
    # G0 = SelfEnergy.G0wrapped(Euv, rtol, kernel.kgrid, param)
    # Σ = SelfEnergy.calcΣ(G0, kernel)

    Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, 1e-7 * param.kF, order, :rpa)
    Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)
    ωgrid = Σ.dlrGrid

    kgrid = Σ.spaceGrid
    kF = kgrid.panel[3]
    kF_label = searchsortedfirst(kgrid.grid, kF)

    ΣR = real(Σ.dynamic)
    ΣI = imag(Σ.dynamic)
    println(Σ.instant[1, 1, :])
    println(ΣR[1, 1, kF_label, :])
    println(ΣI[1, 1, kF_label, :])

    for (n, ω) in enumerate(ΣI[1, 1, kF_label, :])
        println(n, ' ', ωgrid.n[n], ' ', ωgrid.ωn[n], ' ', ω)
    end
    println(SelfEnergy.zfactor(Σ))
    @test isapprox(SelfEnergy.zfactor(Σ), 0.862; rtol = 1e-3)

    G = SelfEnergy.Gwrapped(Σ, param)
    println(G.dynamic[1, 1, kF_label, :])
end
