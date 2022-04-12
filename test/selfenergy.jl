@testset "Self Energy" begin

    @testset "3D Fock" begin
        θ, rs = 0.1, 1.0
        para = Parameter.rydbergUnit(θ, rs, 3)
        println("$(para.μ), $(para.EF)")

        factor = -para.e0^2 * para.kF / π
        #test edge case when k → 0
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(0.0, para) / factor, 2.0, rtol=1e-6)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(1e-6, para) / factor, 2.0, rtol=1e-6)

        #test edge case when k → 0
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para.kF, para) / factor, 1.0, rtol=1e-6)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para.kF + 1e-7, para) / factor, 1.0, rtol=1e-6)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para.kF - 1e-7, para) / factor, 1.0, rtol=1e-6)

        #test edge case when Λs → 0
        para = Parameter.rydbergUnit(θ, rs, 3, Λs=1e-12)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(0.0, para) / factor, 2.0, rtol=1e-6)
    end

    @testset "2D Fock" begin
        θ, rs = 0.1, 1.0
        para = Parameter.rydbergUnit(θ, rs, 2, Λs=0.1)
        println("$(para.μ), $(para.EF)")
        #test edge case when k → 0
        f1 = SelfEnergy.Fock0_ZeroTemp(0.0, para)
        f2 = SelfEnergy.Fock0_ZeroTemp(1e-7, para)
        @test isapprox(f1, f2, rtol=1e-6)

        para1 = Parameter.rydbergUnit(θ, rs, 2)
        para2 = Parameter.rydbergUnit(θ, rs, 2, Λs=1e-12)
        @test isapprox(SelfEnergy.Fock0_ZeroTemp(para1.kF, para1), SelfEnergy.Fock0_ZeroTemp(para2.kF, para2), rtol=1e-6)
    end

    @testset "3D RPA" begin
        # make sure everything works for different unit sets
        θ = 1e-5
        rslist = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        zlist = [0.859, 0.764, 0.700, 0.645, 0.602, 0.568]
        mlist = [0.970, 0.992, 1.016, 1.039, 1.059, 1.078]
        ## our results: zlist = [0.859, 0.764, 0.693, 0.637, 0.591]
        # rslist = [1.0, 2.0]
        # zlist = [0.859, 0.764]
        # mlist = [0.970, 0.992]

        for (ind, rs) in enumerate(rslist)
            param = Parameter.rydbergUnit(θ, rs)
            kF, β = param.kF, param.β
            Euv, rtol = 1000 * param.EF, 1e-11
            # Nk, order = 8, 4
            # maxK, minK = 10kF, 1e-7kF
            # Nk, order = 11, 8
            maxK, minK = 20kF, 1e-8kF
            Nk, order = 16, 12
            # maxK, minK = 10kF, 1e-9kF
            # Euv, rtol = 1000 * param.EF, 1e-11
            # maxK, minK = 20param.kF, 1e-9param.kF
            # Nk, order = 16, 12

            # test G0
            kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [0.0, kF], Nk, minK, order)
            G0 = SelfEnergy.G0wrapped(Euv, rtol, kgrid, param)
            kF_label = searchsortedfirst(kgrid.grid, kF)
            G_tau = GreenFunc.toTau(G0)
            # println(G_tau.timeGrid.grid)
            # println(real(G_tau.dynamic[1, 1, :, end]) .* (-1))
            G_ins = tau2tau(G_tau.dlrGrid, G_tau.dynamic, [β,], G_tau.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
            # println(real(G_ins))
            integrand = real(G_ins) .* kgrid.grid .* kgrid.grid
            density0 = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π^2
            # println("$density0, $(param.n)")
            # @test isapprox(param.n, density0, rtol=3e-5)

            @time Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)
            Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

            Z0 = (SelfEnergy.zfactor(Σ))
            z = zlist[ind]
            # @test isapprox(Z0, z, rtol=3e-3)
            mratio = SelfEnergy.massratio(param, Σ)
            m = mlist[ind]
            # @test isapprox(mratio, m, rtol=3e-3)

            println("θ = $θ,  rs= $rs")
            println("Z-factor = $Z0 ($z), rtol=$(Z0/z-1)")
            println("m*/m = $mratio ($m), rtol=$(mratio/m-1)")

            ## test G_RPA
            G = SelfEnergy.Gwrapped(Σ, param)
            kF_label = searchsortedfirst(G.spaceGrid.grid, kF)
            G_tau = GreenFunc.toTau(G)
            G_ins = tau2tau(G_tau.dlrGrid, G_tau.dynamic, [β,], G_tau.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
            # println(real(G_tau.dynamic[1, 1, kF_label, :]))
            # println(real(G_ins[kF_label]))
            Gk_ins = CompositeGrids.Interp.interp1DGrid(G_ins, G.spaceGrid, kgrid.grid)
            integrand = real(Gk_ins) .* kgrid.grid .* kgrid.grid
            density = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π^2
            # @test isapprox(param.n, density, rtol=2e-4)
            println("density n, from G0, from G_RPA: $(param.n),  $density0 (rtol=$(density0/param.n-1)),  $density (rtol=$(density/param.n-1))")
        end
    end

    @testset "2D RPA" begin
        dim = 2
        θ = 1e-5
        rslist = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0]
        zlist = [0.786, 0.662, 0.519, 0.437, 0.383, 0.344, 0.270, 0.240]
        mlist = [0.981, 1.020, 1.078, 1.117, 1.143, 1.162, 1.196, 1.209]
        # rslist = [0.5, 1.0, 2.0]
        # zlist = [0.786, 0.662, 0.519]
        # mlist = [0.981, 1.020, 1.078]
        for (ind, rs) in enumerate(rslist)
            param = Parameter.rydbergUnit(θ, rs, dim)
            kF, β = param.kF, param.β

            Euv, rtol = 1000 * param.EF, 1e-10
            # set Nk, minK = 8, 1e-7 for β<1e6;  11, 1e-8 for β<1e7
            # Nk, order, minK = 11, 8, 1e-8
            Nk, order = 8, 8
            maxK, minK = 10kF, 1e-7kF
            # Euv, rtol = 100 * param.EF, 1e-12
            # maxK, minK = 20param.kF, 1e-9param.kF
            # Nk, order = 16, 6

            ## test G0
            kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order)
            G0 = SelfEnergy.G0wrapped(Euv, rtol, kgrid, param)
            kF_label = searchsortedfirst(kgrid.grid, kF)
            G_tau = GreenFunc.toTau(G0)
            G_ins = tau2tau(G_tau.dlrGrid, G_tau.dynamic, [β,], G_tau.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
            integrand = real(G_ins) .* kgrid.grid
            density0 = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π
            # @test isapprox(param.n, density0, rtol=3e-5)

            @time Σ = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)
            Σ = SelfEnergy.GreenFunc.toMatFreq(Σ)

            kgrid = Σ.spaceGrid
            kF = kgrid.panel[3]
            Z0 = (SelfEnergy.zfactor(Σ))
            z = zlist[ind]
            # @test isapprox(Z0, z, rtol=3e-3)
            mratio = SelfEnergy.massratio(param, Σ)
            m = mlist[ind]
            # @test isapprox(mratio, m, rtol=3e-3)
            println("θ = $θ,  rs= $rs")
            println("Z-factor = $Z0 ($z), rtol=$(Z0/z-1)")
            println("m*/m = $mratio ($m), rtol=$(mratio/m-1)")

            ## test G_RPA
            G = SelfEnergy.Gwrapped(Σ, param)
            kF_label = searchsortedfirst(kgrid.grid, kF)
            G_tau = GreenFunc.toTau(G)
            G_ins = tau2tau(G_tau.dlrGrid, G_tau.dynamic, [β,], G_tau.timeGrid.grid; axis=4)[1, 1, :, 1] .* (-1)
            # println(real(G_tau.dynamic[1, 1, kF_label, :]))
            # println(real(G_ins[kF_label]))
            # println(real(G_ins))
            integrand = real(G_ins) .* kgrid.grid
            density = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π
            # @test isapprox(param.n, density, rtol=2e-4)
            println("density n, from G0, from G_RPA: $(param.n),  $density0 (rtol=$(density0/param.n-1)),  $density (rtol=$(density/param.n-1))")
        end
    end
end