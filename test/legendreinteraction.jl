@testset "LegendreInteraction" begin

    # set parameters
    param = LegendreInteraction.Parameter.defaultUnit(1/1000.0, 1.0)
    kF = param.kF
    Nk, minK, order, maxK = 16, 1e-8, 6, 10.0
    print(param)

    @testset "Helper function" begin
        # test helper function
        # test with known result: bare coulomb interaction
        H1=LegendreInteraction.helper_function(
            1.0, 1, u->LegendreInteraction.interaction_instant(u,param,:sigma), param;
            Nk=Nk,minK=minK,order=order
        )
        H2=LegendreInteraction.helper_function(
            2.0, 1, u->LegendreInteraction.interaction_instant(u,param,:sigma), param;
            Nk=Nk,minK=minK,order=order
        )
        # println("$(H2-H1) == $(4*π*param.e0^2*log(2))")
        @test isapprox(H2-H1, 4*π*param.e0^2*log(2), rtol = 1e-4)

        helper_grid = CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.1*maxK], [0.0, 2kF], 4, 0.001, 4)
        intgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, helper_grid[end]], [0.0,2kF], 16Nk, 0.01minK, 8order)
        helper = LegendreInteraction.helper_function_grid(helper_grid,intgrid, 1, u->LegendreInteraction.interaction_instant(u,param,:sigma),param)
        helper_analytic = (4*π*param.e0^2) .* log.(helper_grid.grid)
        helper_old = zeros(Float64, helper_grid.size)
        for (yi, y) in enumerate(helper_grid)
            helper_old[yi] = LegendreInteraction.helper_function(y, 1, u->LegendreInteraction.interaction_instant(u,param,:sigma),param)
        end
        # println(helper .- helper[1])
        # println(helper_old .- helper_old[1])
        # println(helper_analytic .- helper_analytic[1])
        rtol = 1e-6
        for (yi, y) in enumerate(helper_grid)
            @test isapprox(helper[yi]-helper[1], helper_analytic[yi]-helper_analytic[1], rtol=rtol)
            @test isapprox(helper_old[yi]-helper_old[1], helper_analytic[yi]-helper_analytic[1], rtol=rtol)
        end
    end

    @testset "DCKernel" begin
        # test DCKernel

        kernel = LegendreInteraction.DCKernel0(param; spin_state=:sigma, Nk=8, rtol=1e-10)
        # kernel = LegendreInteraction.DCKernel(param, 100*param.EF, 1e-8, 5, 10*param.kF, 1e-7*param.kF, 4, :rpa,0,:sigma)
        # kernel2 = LegendreInteraction.DCKernel(kernel, 500.0)
        kF_label = searchsortedfirst(kernel.kgrid.grid, kernel.param.kF)
        qF_label = searchsortedfirst(kernel.qgrids[kF_label].grid, kernel.param.kF)
        println(kernel.kernel[kF_label,qF_label,:])
    end
end
