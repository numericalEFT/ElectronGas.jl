@testset "Self Energy" begin

    @testset "3D Fock" begin
        θ, rs = 0.1, 1.0
        para = Parameter.rydbergUnit(θ, rs, 3)
        factor = -para.e0^2 * para.kF / π
        #test edge case when k → 0
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(0.0, para) / factor, 2.0, rtol = 1e-6)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(1e-6, para) / factor, 2.0, rtol = 1e-6)

        #test edge case when k → 0
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para.kF, para) / factor, 1.0, rtol = 1e-6)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para.kF + 1e-7, para) / factor, 1.0, rtol = 1e-6)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para.kF - 1e-7, para) / factor, 1.0, rtol = 1e-6)

        #test edge case when Λs → 0
        para = Parameter.rydbergUnit(θ, rs, 3, Λs = 1e-12)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(0.0, para) / factor, 2.0, rtol = 1e-6)
    end

    @testset "2D Fock" begin
        θ, rs = 0.1, 1.0
        para = Parameter.rydbergUnit(θ, rs, 2, Λs = 0.1)
        #test edge case when k → 0
        f1 = SelfEnergy.Fock0_ZeroTemp(0.0, para)
        f2 = SelfEnergy.Fock0_ZeroTemp(1e-7, para)
        @test isapprox(f1, f2, rtol = 1e-6)

        para1 = Parameter.rydbergUnit(θ, rs, 2)
        para2 = Parameter.rydbergUnit(θ, rs, 2, Λs = 1e-12)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para1.kF, para1), SelfEnergy.Fock0_ZeroTemp(para2.kF, para2), rtol = 1e-6)
    end

    @testset "3D RPA" begin
        # make sure everything works for different unit sets
        θ, rs = 1e-3, 1.0
        param = Parameter.defaultUnit(θ, rs)
        Euv, rtol = 100 * param.EF, 1e-10
        Nk, order = 8, 4

        @time Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, 1e-7 * param.kF, order, :rpa)
        Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

        Z0 = (SelfEnergy.zfactor(Σ))
        @test isapprox(Z0, 0.862, rtol = 1e-3)
        println("z-factor = $Z0")
    end

    @testset "2D RPA" begin
        dim = 2
        θ, rs = 1e-3, 1.0
        # Λs = 1e-6
        # param = Parameter.defaultUnit(θ, rs, dim; Λs = Λs)
        param = Parameter.rydbergUnit(θ, rs, dim)

        Euv, rtol = 100 * param.EF, 1e-10
        # set Nk, minK = 8, 1e-7 for β<1e6;  11, 1e-8 for β<1e7
        # Nk, order, minK = 11, 4, 1e-8
        Nk, order, minK = 8, 8, 1e-7

        @time Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, 10 * param.kF, minK * param.kF, order, :rpa)
        Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

        kgrid = Σ.spaceGrid
        kF = kgrid.panel[3]
        Z0 = (SelfEnergy.zfactor(Σ))
        @test isapprox(Z0, 0.662, rtol = 5e-3)

        println("θ = $θ,  rs= $rs")
        println("Z-factor = $Z0")
        # G = SelfEnergy.Gwrapped(Σ, param)
        # println(G.dynamic[1, 1, kF_label, :])
    end
end
