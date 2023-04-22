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

    @testset "RPA based on the user-defined kgrid" begin
        θ, rs = 0.01, 4.0
        para = Parameter.rydbergUnit(θ, rs, 3, Λs=1e-5)
        sigma1 = SelfEnergy.G0W0(para)
        kFidx = locate(sigma1[1].mesh[2], para.kF)
        sigma2 = SelfEnergy.G0W0(para, [sigma1[1].mesh[2][kFidx],])
        @test isapprox(sigma1[1][1, kFidx], sigma2[1][1, 1], rtol=1e-6)
        @test isapprox(sigma1[2][1, kFidx], sigma2[2][1, 1], rtol=1e-6)

        para = Parameter.rydbergUnit(θ, rs, 2, Λs=1e-5)
        sigma1 = SelfEnergy.G0W0(para)
        kFidx = locate(sigma1[1].mesh[2], para.kF)
        sigma2 = SelfEnergy.G0W0(para, [sigma1[1].mesh[2][kFidx],])
        @test isapprox(sigma1[1][1, kFidx], sigma2[1][1, 1], rtol=1e-6)
        @test isapprox(sigma1[2][1, kFidx], sigma2[2][1, 1], rtol=1e-6)
    end

    @testset "3D RPA" begin
        # make sure everything works for different unit sets
        θ = 0.001
        # rslist = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        # zlist = [0.859, 0.764, 0.700, 0.645, 0.602, 0.568]
        # mlist = [0.970, 0.992, 1.016, 1.039, 1.059, 1.078]
        ## our results: zlist = [0.859, 0.764, 0.693, 0.637, 0.591]
        # rslist = [5.0,]
        # zlist = [0.5913,]
        # mlist = [1.059,]
        # rslist = [1.0, 2.0]
        # zlist = [0.859, 0.764]
        # mlist = [0.970, 0.992]
        rslist = [2.0,]
        zlist = [0.764,]
        mlist = [0.992,]

        for (ind, rs) in enumerate(rslist)
            param = Parameter.rydbergUnit(θ, rs)
            kF, β = param.kF, param.β
            Euv, rtol = 1000 * param.EF, 1e-11
            # Nk, order = 8, 4
            # maxK, minK = 10kF, 1e-7kF
            # Nk, order = 11, 8
            maxK, minK = 20kF, 1e-8kF
            Nk, order = 12, 8
            # maxK, minK = 10kF, 1e-9kF
            # Euv, rtol = 1000 * param.EF, 1e-11
            # maxK, minK = 20param.kF, 1e-9param.kF
            # Nk, order = 16, 12
            kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [0.0, kF], Nk, minK, order)

            # test G0
            G0 = SelfEnergy.G0wrapped(Euv, rtol, kgrid, param)
            # kF_label = searchsortedfirst(kgrid.grid, kF)
            kF_label = locate(G0.mesh[2], kF)

            G0_dlr = G0 |> to_dlr
            G0_tau = G0_dlr |> to_imtime

            G0_ins = dlr_to_imtime(G0_dlr, [β,]) * (-1)
            integrand = real(G0_ins[1, :]) .* kgrid.grid .* kgrid.grid
            density0 = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π^2
            @test isapprox(param.n, density0, rtol=3e-5)

            @time Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)

            Z0 = (SelfEnergy.zfactor(param, Σ))[1]
            z = zlist[ind]
            @test isapprox(Z0, z, rtol=3e-3)
            mratio = SelfEnergy.massratio(param, Σ, Σ_ins)[1]
            mratio1 = SelfEnergy.massratio(param, Σ, Σ_ins, 1e-5)[1]
            m = mlist[ind]
            @test isapprox(mratio, m, rtol=3e-3)
            @test isapprox(mratio, mratio1, rtol=1e-4)

            println("θ = $θ,  rs= $rs")
            println("Z-factor = $Z0 ($z), rtol=$(Z0/z-1)")
            println("m*/m = $mratio ($m), rtol=$(mratio/m-1)")

            ## test G_RPA
            G = SelfEnergy.Gwrapped(Σ, Σ_ins, param)
            G_dlr = G |> to_dlr
            G_tau = G_dlr |> to_imtime

            G_ins = dlr_to_imtime(G_dlr, [β,]) * (-1)
            integrand = real(G_ins[1, :]) .* kgrid.grid .* kgrid.grid
            density = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π^2
            @test isapprox(param.n, density, rtol=5e-3)
            println("density n, from G0, from G_RPA: $(param.n),  $density0 (rtol=$(density0/param.n-1)),  $density (rtol=$(density/param.n-1))")
        end
    end

    @testset "2D RPA" begin
        dim = 2
        θ = 1e-5
        # rslist = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0]
        # zlist = [0.786, 0.662, 0.519, 0.437, 0.383, 0.344, 0.270, 0.240]
        # mlist = [0.981, 1.020, 1.078, 1.117, 1.143, 1.162, 1.196, 1.209]
        rslist = [1.0,]
        zlist = [0.662,]
        mlist = [1.02,]
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
            kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order)

            ## test G0
            G0 = SelfEnergy.G0wrapped(Euv, rtol, kgrid, param)
            kF_label = locate(G0.mesh[2], kF)

            G0_dlr = G0 |> to_dlr
            G0_tau = G0_dlr |> to_imtime

            G0_ins = dlr_to_imtime(G0_dlr, [β,]) * (-1)
            integrand = real(G0_ins[1, :]) .* kgrid.grid
            density0 = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π
            @test isapprox(param.n, density0, rtol=3e-5)

            @time Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, :rpa)

            Z0 = (SelfEnergy.zfactor(param, Σ))[1]
            z = zlist[ind]
            @test isapprox(Z0, z, rtol=3e-3)
            mratio = SelfEnergy.massratio(param, Σ, Σ_ins)[1]
            m = mlist[ind]
            @test isapprox(mratio, m, rtol=3e-3)
            println("θ = $θ,  rs= $rs")
            println("Z-factor = $Z0 ($z), rtol=$(Z0/z-1)")
            println("m*/m = $mratio ($m), rtol=$(mratio/m-1)")

            ## test G_RPA
            G = SelfEnergy.Gwrapped(Σ, Σ_ins, param)
            G_dlr = G |> to_dlr
            G_tau = G_dlr |> to_imtime

            G_ins = dlr_to_imtime(G_dlr, [β,]) * (-1)
            integrand = real(G_ins[1, :]) .* kgrid.grid
            density = CompositeGrids.Interp.integrate1D(integrand, kgrid) / π
            @test isapprox(param.n, density, atol=3e-4)
            println("density n, from G0, from G_RPA: $(param.n),  $density0 (rtol=$(density0/param.n-1)),  $density (rtol=$(density/param.n-1))")
        end
    end
end