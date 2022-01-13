@testset "LegendreInteraction" begin

    @testset "DCKernel" begin
        param = LegendreInteraction.Parameter.defaultUnit(1000.0, 1.0)
        print(param)
        kernel = LegendreInteraction.DCKernel0(param; spin_state=:sigma, Nk=8, rtol=1e-10)
        # kernel = LegendreInteraction.DCKernel(param, 100*param.EF, 1e-8, 5, 10*param.kF, 1e-7*param.kF, 4, :rpa,0,:sigma)
        # kernel2 = LegendreInteraction.DCKernel(kernel, 500.0)
        kF_label = searchsortedfirst(kernel.kgrid.grid, kernel.param.kF)
        qF_label = searchsortedfirst(kernel.qgrids[kF_label].grid, kernel.param.kF)
        println(kernel.kernel[kF_label,qF_label,:])
        Nk, mink, order = 16, 1e-8, 6
        H1=LegendreInteraction.helper_function(
            1.0, 1, u->LegendreInteraction.interaction_instant(u,param,:sigma), param;
            Nk=Nk,mink=mink,order=order
        )
        H2=LegendreInteraction.helper_function(
            2.0, 1, u->LegendreInteraction.interaction_instant(u,param,:sigma), param;
            Nk=Nk,mink=mink,order=order
        )
        println("$(H2-H1) == $(4*π*param.e0^2*log(2))")
    end
end
