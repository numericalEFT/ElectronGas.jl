# using Gaston

@testset "Polarization" begin
    dim = 2
    beta, rs = 1e2, 1.0
    param = Parameter.defaultUnit(1 / beta, rs, dim)
    @testset "Single Point Polarization" begin
        # for low temp, Π from zero temp and finite temp should be close
        testq = [-1.0, 0.0, 1e-160, 1e-8, 0.5, 1.0, 2.0, 10.0]
        testn = [0, 1, 100]

        for q in testq
            for n in testn
                PZ = Polarization.Polarization0_ZeroTemp(q, n, param)
                PF = Polarization.Polarization0_FiniteTemp(q, n, param)
                # @test abs((PZ - PF) / (PZ + 1e-6)) < 1e-6
                # println("q=$q, n=$n")
                # println(Polarization.Polarization0_ZeroTemp(q, n, param), "\n")
                # println(Polarization.Polarization0_FiniteTemp(q, n, param))
            end
        end
    end

    @testset "Polarization vs. q" begin
        qgrid = [0.0, 1e-160, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]
        polar = similar(qgrid)
        n = 2
        println("n = $n,  ω_n = $(2π*n/param.β)")
        for (qi, q) in enumerate(qgrid)
            polar[qi] = Polarization.Polarization0_ZeroTemp(q, n, param)
            # println("q=$q, Π0=$(polar[qi])")
        end
    end

    @testset "Polarization vs. n" begin
        q = 1.0
        ngrid = [n for n in -5:20]
        polar = zeros(Float64, size(ngrid))

        println("q = $q")
        for (i, n) in enumerate(ngrid)
            polar[i] = Polarization.Polarization0_ZeroTemp(q, n, param)
            # println("n=$n, Π0=$(polar[i]) \n")
        end
    end

    @testset "Polarization: ZeroTemp vs. FiniteTemp" begin
        dim = 3
        beta, rs = 1e8, 1.0
        param = Parameter.defaultUnit(1 / beta, rs, dim)

        qgrid = [0.0, 1e-16, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]
        ngrid = [n for n in -5:20]
        for (qi, q) in enumerate(qgrid)
            for (ni, n) in enumerate(ngrid)
                # println("q=$q, n=$n")
                @test isapprox(
                    Polarization.Polarization0_ZeroTemp(q, n, param),
                    Polarization.Polarization0_FiniteTemp(q, n, param),
                    rtol = 1e-6
                )
            end
        end

        # dim = 2
        # beta, rs = 1e8, 1.0
        # param = Parameter.defaultUnit(1 / beta, rs, dim)

        # qgrid = [0.0, 1e-160, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]
        # ngrid = [n for n in -5:20]
        # for (qi, q) in enumerate(qgrid)
        #     for (ni, n) in enumerate(ngrid)
        #         @test isapprox(
        #             Polarization.Polarization0_ZeroTemp(q, n, param),
        #             Polarization.Polarization0_FiniteTemp(q, n, param),
        #             rtol = 1e-4
        #         )
        #     end
        # end

    end

    @testset "Polarization: wrapped" begin
        dim = 3
        beta, rs = 1e3, 1.0
        param = Parameter.defaultUnit(1 / beta, rs, dim)

        qgrid = [0.0, 1e-16, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]

        PZ = Polarization.Polarization0wrapped(100param.EF, 1e-10, qgrid, param, pifunc = Polarization.Polarization0_ZeroTemp)
        PF = Polarization.Polarization0wrapped(100param.EF, 1e-10, qgrid, param, pifunc = Polarization.Polarization0_FiniteTemp)

    end

    # while true
    #     print("Please enter a whole number between 1 and 5: ")
    #     input = readline(stdin)
    #     value = tryparse(Int, input)
    #     if value !== nothing && 1 <= value <= 5
    #         println("You entered $(input)")
    #         break
    #     else
    #         @warn "Enter a whole number between 1 and 5"
    #     end
    # end
    # plt = plot(qgrid, polar, ytics = -0.08:0.01, ls = 1, Axes(grid = :on, key = "left"))
    # display(plt)
    # save(term = "png", output = "polar_2d.png",
    # saveopts = "font 'Consolas,10' size 1280,900 lw 1")
end
