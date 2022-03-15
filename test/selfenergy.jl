@testset "Self Energy" begin

    @testset "3D Fock" begin
        θ, rs = 0.1, 1.0
        para = Parameter.rydbergUnit(θ, rs, 3)
        println("$(para.μ), $(para.EF)")

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
        println("$(para.μ), $(para.EF)")
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
        θ = 1e-4
        # rslist = [1.0, 2.0, 3.0, 4.0, 5.0]
        # zlist = [0.859, 0.764, 0.700, 0.645, 0.602]
        # mlist = [0.970, 0.992, 1.016, 1.039, 1.059]
        ## our results: zlist = [0.859, 0.764, 0.693, 0.637, 0.591]

        rslist = [1.0, 2.0]
        zlist = [0.859, 0.764]
        mlist = [0.970, 0.992]

        for (ind, rs) in enumerate(rslist)
            param = Parameter.defaultUnit(θ, rs)
            Euv, rtol = 100 * param.EF, 1e-10
            # Nk, order = 8, 4
            Nk, order = 11, 8
            maxK, minK = 10param.kF, 1e-8 * param.kF

            @time Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)
            Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

            Z0 = (SelfEnergy.zfactor(Σ))
            z = zlist[ind]
            @test isapprox(Z0, z, atol = 3e-3)
            println("θ = $θ,  rs= $rs")
            println("Z-factor = $Z0 ($z)")

            mratio = SelfEnergy.massratio(param, Σ)
            m = mlist[ind]
            @test isapprox(mratio, m, atol = 3e-3)
            println("m*/m = $mratio ($m)")
        end
    end

    @testset "2D RPA" begin
        dim = 2
        θ = 1e-4
        # rslist = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0]
        # zlist = [0.786, 0.662, 0.519, 0.437, 0.383, 0.344, 0.270, 0.240]
        # mlist = [0.981, 1.020, 1.078, 1.117, 1.143, 1.162, 1.196, 1.209]
        rslist = [0.5, 1.0]
        zlist = [0.786, 0.662]
        mlist = [0.981, 1.020]
        for (ind, rs) in enumerate(rslist)
            param = Parameter.rydbergUnit(θ, rs, dim)

            Euv, rtol = 100 * param.EF, 1e-10
            # set Nk, minK = 8, 1e-7 for β<1e6;  11, 1e-8 for β<1e7
            # Nk, order, minK = 11, 8, 1e-8
            Nk, order = 8, 8
            maxK, minK = 10param.kF, 1e-7param.kF
            # Euv, rtol = 100 * param.EF, 1e-12
            # maxK, minK = 20param.kF, 1e-9param.kF
            # Nk, order = 16, 6

            @time Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)
            Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

            kgrid = Σ.spaceGrid
            kF = kgrid.panel[3]
            Z0 = (SelfEnergy.zfactor(Σ))
            z = zlist[ind]
            @test isapprox(Z0, z, atol = 3e-3)
            println("θ = $θ,  rs= $rs")
            println("Z-factor = $Z0 ($z)")

            mratio = SelfEnergy.massratio(param, Σ)
            m = mlist[ind]
            @test isapprox(mratio, m, atol = 3e-3)
            println("m*/m = $mratio ($m)")
            # G = SelfEnergy.Gwrapped(Σ, param)
            # println(G.dynamic[1, 1, kF_label, :])
        end
    end
end
