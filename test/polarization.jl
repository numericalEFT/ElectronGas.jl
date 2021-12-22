@testset "Polarization" begin

    @testset "Single Point Polarization" begin
	      beta = 1e8
        # for low temp, Π from zero temp and finite temp should be close
        param = Parameter.defaultUnit(beta, 1.0)
        testq = [-1.0, 0.0, 1e-160, 1e-8, 0.5, 1.0, 2.0, 10.0]
        testn = [0, 1, 100]

        for q in testq
            for n in testn
                PZ = Polarization.Polarization0_ZeroTemp(q, n, param)
                PF = Polarization.Polarization0_FiniteTemp(q, n, param)
                @test abs((PZ-PF)/(PZ+1e-6)) < 1e-6
                # println("q=$q, n=$n")
                # println(Polarization.Polarization0_ZeroTemp(q, n, param))
                # println(Polarization.Polarization0_FiniteTemp(q, n, param))
            end
        end

    end
end
