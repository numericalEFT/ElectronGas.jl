@testset "Self Energy" begin
    T= 1/10000.0
    rs = 10.0
    @testset "default unit" begin
        param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(T, rs)
        Euv, rtol = 100*param.EF, 1e-10
        Nk, order = 12, 5

        # kernel = SelfEnergy.LegendreInteraction.DCKernel(param, Euv, rtol, 10, 10*param.kF, 1e-7*param.kF, 5, :rpa,0,:sigma)
        # G0 = SelfEnergy.G0wrapped(Euv, rtol, kernel.kgrid, param)
        # Σ = SelfEnergy.calcΣ(G0, kernel)

        Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10*param.kF, 1e-8*param.kF, order, :rpa)
        Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

        kgrid = Σ.spaceGrid
        kF = kgrid.panel[3]
        kF_label = searchsortedfirst(kgrid.grid, kF)

        ΣR = real(Σ.dynamic)
        ΣI = imag(Σ.dynamic)
        println(Σ.instant[1,1,:])
        println(ΣR[1,1,kF_label,:])
        println(ΣI[1,1,kF_label,:])
        Z0 = (SelfEnergy.zfactor(Σ))
        #@test isapprox(Z0, 0.862, rtol=1e-3)
        println("z-factor = $Z0")
        G = SelfEnergy.Gwrapped(Σ,param)
        println(G.dynamic[1,1,kF_label,:])
    end
    @testset "rydberg unit" begin
        # make sure everything works for different unit sets
        param = SelfEnergy.LegendreInteraction.Parameter.rydbergUnit(T, rs)
        Euv, rtol = 100*param.EF, 1e-10
        Nk, order = 8, 4

        # kernel = SelfEnergy.LegendreInteraction.DCKernel(param, Euv, rtol, 10, 10*param.kF, 1e-7*param.kF, 5, :rpa,0,:sigma)
        # G0 = SelfEnergy.G0wrapped(Euv, rtol, kernel.kgrid, param)
        # Σ = SelfEnergy.calcΣ(G0, kernel)

        Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10*param.kF, 1e-7*param.kF, order, :rpa)
        Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)
        Z0 = (SelfEnergy.zfactor(Σ))
        #@test isapprox(Z0, 0.862, rtol=1e-3)
        println("z-factor = $Z0")
    end
end
