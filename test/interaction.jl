@testset "Interaction" begin
    beta = 400
    rs = 4.0
    param = Parameter.defaultUnit(1 / beta, rs)
    @testset "bubble dyson" begin
        @test Interaction.coulombinv(0.0, Parameter.derive(param, gs=0.0, ga=0.0, Λs=0.0, Λa=0.0)) == (Inf, Inf)
        @test Interaction.coulomb(0.0, Parameter.derive(param, gs=0.0, ga=0.0, Λs=0.0, Λa=0.0)) == (0.0, 0.0)
        @test Interaction.coulombinv(0.0, Parameter.derive(param, ga=1.0, Λs=0.0, Λa=0.0)) == (0.0, 0.0)
        @test Interaction.coulomb(0.0, Parameter.derive(param, ga=1.0, Λs=0.0, Λa=0.0)) == (Inf, Inf)
        @test Interaction.coulombinv(1.0, Parameter.derive(param, gs=0.0, ga=0.0, Λs=0.0, Λa=0.0)) == (Inf, Inf)
        @test Interaction.coulomb(1.0, Parameter.derive(param, gs=0.0, ga=0.0, Λs=0.0, Λa=0.0)) == (0.0, 0.0)

        @test Interaction.bubbledyson(Inf, 1.0, 1.0) == 0.5
        @test Interaction.bubbledyson(Inf, 0.0, 1.0) == 0.0
        @test Interaction.bubbledyson(0.0, 1.0, 1.0) == -Inf
        @test Interaction.bubbledyson(0.0, 0.0, 1.0) == -Inf

    end

    @testset "RPA and KO" begin
        testq = [-1.0, 0.0, 1e-160, 1e-8, 0.5, 1.0, 2.0, 10.0]
        val = Interaction.RPA(1.0, 1, param)[1] / Interaction.coulomb(1.0, param)[1]
        @test isapprox(Interaction.RPA(1.0, 1, param; regular=true)[1], val, rtol=1e-10)
        println(Interaction.KO(1.0, 1, param))
        for q in testq
            println(Interaction.RPA(q, 1, param; regular=true))
        end
        println(Interaction.KO(1.0, 1, param; regular=true))
        # println(Interaction.RPA(1.0, 1, param; pifunc = Interaction.Polarization.Polarization0_FiniteTemp))
        # println(Interaction.KO(1.0, 1, param; pifunc = Interaction.Polarization.Polarization0_FiniteTemp))

        RPAs, RPAa = Interaction.RPAwrapped(100 * param.EF, 1e-8, [1e-8, 0.5, 1.0, 2.0, 10.0], param)
        # println(RPAs.dynamic)
        # println(RPAs.instant)
        KOs, KOa = Interaction.KOwrapped(100 * param.EF, 1e-8, [1e-8, 0.5, 1.0, 2.0, 10.0], param)
        # println(KOs.dynamic)
        # println(KOs.instant)
    end
end
