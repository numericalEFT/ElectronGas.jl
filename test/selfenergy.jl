@testset "Self Energy" begin
    @testset "default unit" begin
        dim = 2
        θ, rs = 1e-3, 1.0
        param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(θ, rs, dim)

        Euv, rtol = 100 * param.EF, 1e-10
        # set Nk, minK = 8, 1e-7 for β<1e6;  11, 1e-8 for β<1e7
        Nk, order, minK = 11, 4, 1e-8

        Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, minK * param.kF, order, :rpa)
        Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

        kgrid = Σ.spaceGrid
        kF = kgrid.panel[3]
        kF_label = searchsortedfirst(kgrid.grid, kF)

        ΣR = real(Σ.dynamic)
        ΣI = imag(Σ.dynamic)
        println(Σ.instant[1, 1, :])
        println(ΣR[1, 1, kF_label, :])
        println(ΣI[1, 1, kF_label, :])
        Z0 = (SelfEnergy.zfactor(Σ))
        # @test isapprox(Z0, 0.862, rtol = 5e-3)
        println("Z-factor = $Z0")
        G = SelfEnergy.Gwrapped(Σ, param)
        println(G.dynamic[1, 1, kF_label, :])
    end
    # @testset "default unit" begin
    #     # make sure everything works for different unit sets
    #     param = SelfEnergy.LegendreInteraction.Parameter.rydbergUnit(1 / 1000.0, 1.0)
    #     Euv, rtol = 100 * param.EF, 1e-10
    #     Nk, order = 8, 4

    #     # kernel = SelfEnergy.LegendreInteraction.DCKernel(param, Euv, rtol, 10, 10*param.kF, 1e-7*param.kF, 5, :rpa,0,:sigma)
    #     # G0 = SelfEnergy.G0wrapped(Euv, rtol, kernel.kgrid, param)
    #     # Σ = SelfEnergy.calcΣ(G0, kernel)

    #     Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, 1e-7 * param.kF, order, :rpa)
    #     Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)
    #     Z0 = (SelfEnergy.zfactor(Σ))
    #     @test isapprox(Z0, 0.862, rtol = 1e-3)
    #     println("z-factor = $Z0")
    # end
end
