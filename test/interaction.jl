@testset "Interaction" begin
    beta = 1e4
    rs = 2.0
    param = Parameter.defaultUnit(1 / beta, rs)
    @testset "bubble dyson" begin
        @test Interaction.coulombinv(0.0, Parameter.derive(param, e0 = 0.0, espin = 0.0, Λs = 0.0, Λa = 0.0)) == (Inf, Inf)
        @test Interaction.coulomb(0.0, Parameter.derive(param, e0 = 0.0, espin = 0.0, Λs = 0.0, Λa = 0.0)) == (0.0, 0.0)
        @test Interaction.coulombinv(0.0, Parameter.derive(param, e0 = 1.0, espin = 1.0, Λs = 0.0, Λa = 0.0)) == (0.0, 0.0)
        @test Interaction.coulomb(0.0, Parameter.derive(param, e0 = 1.0, espin = 1.0, Λs = 0.0, Λa = 0.0)) == (Inf, Inf)
        @test Interaction.coulombinv(1.0, Parameter.derive(param, e0 = 0.0, espin = 0.0, Λs = 0.0, Λa = 0.0)) == (Inf, Inf)
        @test Interaction.coulomb(1.0, Parameter.derive(param, e0 = 0.0, espin = 0.0, Λs = 0.0, Λa = 0.0)) == (0.0, 0.0)

        @test Interaction.bubbledyson(Inf, 1.0, 1.0) == 0.5
        @test Interaction.bubbledyson(Inf, 0.0, 1.0) == 0.0
        @test Interaction.bubbledyson(0.0, 1.0, 1.0) == -Inf
        @test Interaction.bubbledyson(0.0, 0.0, 1.0) == -Inf

    end

    @testset "RPA and KO" begin
        println(Interaction.RPA(1.0, 1, param))
        println(Interaction.KO(1.0, 1, param))
        println(Interaction.RPA(1.0, 1, param; regular = true))
        println(Interaction.KO(1.0, 1, param; regular = true))
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
