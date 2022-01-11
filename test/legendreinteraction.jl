@testset "LegendreInteraction" begin

    @testset "DCKernel" begin
        param = LegendreInteraction.Parameter.defaultUnit(1000.0, 1.0)

        kernel = LegendreInteraction.DCKernel(param, 100*param.EF, 1e-8, 5, 10*param.kF, 1e-7*param.kF, 4, :rpa,0,:sigma)
        # kernel2 = LegendreInteraction.DCKernel(kernel, 500.0)
        kF_label = searchsortedfirst(kernel.kgrid.grid, kernel.param.kF)
        qF_label = searchsortedfirst(kernel.qgrids[kF_label].grid, kernel.param.kF)
        println(kernel.kernel[kF_label,qF_label,:])
    end
end