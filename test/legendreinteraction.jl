@testset "LegendreInteraction" begin

    # set parameters
    param = LegendreInteraction.Parameter.defaultUnit(1 / 1000.0, 1.0)
    kF = param.kF
    Nk, minK, order, maxK = 16, 1e-8, 6, 10.0
    # print(param)

    @testset "Helper function" begin
        # test helper function
        # test with known result: bare coulomb interaction
        H1 = LegendreInteraction.helper_function(
            1.0, 1, u -> LegendreInteraction.interaction_instant(u, param, :sigma), param;
            Nk = Nk, minK = minK, order = order
        )
        H2 = LegendreInteraction.helper_function(
            2.0, 1, u -> LegendreInteraction.interaction_instant(u, param, :sigma), param;
            Nk = Nk, minK = minK, order = order
        )
        # println("$(H2-H1) == $(4*π*param.e0^2*log(2))")
        @test isapprox(H2 - H1, 4 * π * param.e0^2 * log(2), rtol = 1e-4)

        helper_grid = CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.1 * maxK], [0.0, 2kF], 4, 0.001, 4)
        intgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, helper_grid[end]], [0.0, 2kF], 2Nk, 0.01minK, 2order)
        helper = LegendreInteraction.helper_function_grid(helper_grid, intgrid, 1, u -> LegendreInteraction.interaction_instant(u, param, :sigma), param)
        helper_analytic = (4 * π * param.e0^2) .* log.(helper_grid.grid)
        helper_old = zeros(Float64, helper_grid.size)
        for (yi, y) in enumerate(helper_grid)
            helper_old[yi] = LegendreInteraction.helper_function(y, 1, u -> LegendreInteraction.interaction_instant(u, param, :sigma), param)
        end
        # println(helper .- helper[1])
        # println(helper_old .- helper_old[1])
        # println(helper_analytic .- helper_analytic[1])
        rtol = 1e-6
        for (yi, y) in enumerate(helper_grid)
            @test isapprox(helper[yi] - helper[1], helper_analytic[yi] - helper_analytic[1], rtol = rtol)
            @test isapprox(helper_old[yi] - helper_old[1], helper_analytic[yi] - helper_analytic[1], rtol = rtol)
        end
    end

    @testset "DCKernel_2d" begin
        # set parameters
        param = LegendreInteraction.Parameter.defaultUnit(1e-4, 1.0, dim = 2)
        kF = param.kF
        println("kF = $kF")
        # Euv, rtol = 100 * param.EF, 1e-10
        # Nk, minK, order, maxK = 8, 1e-7kF, 8, 10kF

        Euv, rtol = 100 * param.EF, 1e-12
        maxK, minK = 20param.kF, 1e-9param.kF
        Nk, order = 16, 6
        int_type = :rpa

        W = LegendreInteraction.DCKernel_2d(param;
            Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order, int_type = int_type, spin_state = :sigma)
        fdlr = ElectronGas.Lehmann.DLRGrid(Euv, param.β, rtol, true, :pha)
        bdlr = W.dlrGrid
        kgrid = W.kgrid
        qgrids = W.qgrids

        kernel_bare = W.kernel_bare
        kernel_freq = W.kernel
        kernel = real(ElectronGas.Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis = 3))
        kF_label = searchsortedfirst(kgrid.grid, kF)
        qF_label = searchsortedfirst(qgrids[kF_label].grid, kF)

        println("static kernel at (kF, kF): $(kgrid.grid[kF_label]), $(qgrids[kF_label].grid[qF_label])")
        println("$(kernel_bare[kF_label, qF_label])")
        # println(kernel_bare[kF_label, :])
        @test isapprox(kernel_bare[kF_label, qF_label], 284, rtol = 0.04)

        println("dynamic kernel at (kF, kF):")
        println(view(kernel, kF_label, qF_label, :))
    end

    @testset "Test case: r_s=4, β=400, 3D" begin
        param = Parameter.defaultUnit(1 / 400.0, 4.0) # 3D UEG
        Euv, rtol = 100 * param.EF, 1e-12
        maxK, minK = 20param.kF, 1e-9param.kF
        Nk, order = 16, 6
        int_type = :rpa

        #--- prepare kernel ---
        W = LegendreInteraction.DCKernel0(param;
            Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order,
            int_type = int_type)

        fdlr = ElectronGas.Lehmann.DLRGrid(Euv, param.β, rtol, true, :pha)
        bdlr = W.dlrGrid
        kgrid = W.kgrid
        qgrids = W.qgrids

        kernel_bare = W.kernel_bare
        kernel_freq = W.kernel
        kernel = real(ElectronGas.Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis = 3))
        kF_label = searchsortedfirst(kgrid.grid, param.kF)
        qF_label = searchsortedfirst(qgrids[kF_label].grid, param.kF)

        k, p = kgrid.grid[kF_label], qgrids[kF_label].grid[qF_label]
        w0 = 4π * param.e0^2 * log((k + p) / abs(k - p)) / k / p

        println("static kernel at (kF, kF): $k, $p")
        println("$(kernel_bare[kF_label, qF_label]),  Analytic: $w0")
        @test isapprox(kernel_bare[kF_label, qF_label], w0, rtol = 1e-6)
        # @test isapprox(kernel_bare[kF_label, qF_label], 1314.5721928, rtol = 1e-6)
        # println(kernel_bare[kF_label, :])

        println("dynamic kernel at (kF, kF):")
        println(view(kernel, kF_label, qF_label, :))

    end
end
