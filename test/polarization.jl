using Gaston

@testset "Polarization" begin

    @testset "Single Point Polarization" begin
        dim = 2
        beta, rs = 1e8, 1.0
        # for low temp, Π from zero temp and finite temp should be close
        param = Parameter.defaultUnit(1 / beta, rs, dim)
        testq = [-1.0, 0.0, 1e-160, 1e-8, 0.5, 1.0, 2.0, 10.0]
        testn = [0, 1, 100]

        for q in testq
            for n in testn
                PZ = Polarization.Polarization0_ZeroTemp(q, n, param)
                PF = Polarization.Polarization0_FiniteTemp(q, n, param)
                # @test abs((PZ - PF) / (PZ + 1e-6)) < 1e-6
                println("q=$q, n=$n")
                println(Polarization.Polarization0_ZeroTemp(q, n, param))
                println(Polarization.Polarization0_FiniteTemp(q, n, param))
            end
        end

        # qgrid = [0.0, 1e-160, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.15, 0.3, 0.5, 1.0]
        # polar = similar(qgrid)

        # n = Int(5e6)
        # println("ω_n = $(2π*n/param.β)")
        # for (qi, q) in enumerate(qgrid)
        #     polar[qi] = Polarization.Polarization0_ZeroTemp(q, n, param)
        #     println("q=$q, Π0=$(polar[qi]) \n")
        # end
        # plt = plot(qgrid, polar, ytics = -0.08:0.01, ls = 1, Axes(grid = :on, key = "left"))
        # display(plt)

        while true
            print("Please enter a whole number between 1 and 5: ")
            input = readline(stdin)
            value = tryparse(Int, input)
            if value !== nothing && 1 <= value <= 5
                println("You entered $(input)")
                break
            else
                @warn "Enter a whole number between 1 and 5"
            end
        end
        # save(term = "png", output = "polar_2d.png",
        # saveopts = "font 'Consolas,10' size 1280,900 lw 1")

    end
end
