@testset "Interaction" begin

    @testset "RPA and KO" begin
        beta = 1e4
        rs = 2.0
        param = Interaction.Parameter.defaultUnit(beta,rs)
        println(Interaction.RPA(0.01, 1, param))
        println(Interaction.KO(0.01, 1, param))
        println(Interaction.RPA(1.0, 1, param, Interaction.Polarization.Polarization0_FiniteTemp))
        println(Interaction.KO(1.0, 1, param, Interaction.Polarization.Polarization0_FiniteTemp))

        RPA, = Interaction.RPAwrapped(100*param.EF,1e-8,[1e-8,0.5,1.0,2.0,10.0],param)
        println(RPA.dynamic)
        println(RPA.instant)
        KOs, KOa = Interaction.KOwrapped(100*param.EF,1e-8,[1e-8,0.5,1.0,2.0,10.0],param)
        println(KOs.dynamic)
        println(KOs.instant)
    end
end