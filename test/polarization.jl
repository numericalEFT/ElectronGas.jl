# using Gaston

@testset "Polarization" begin
    beta, rs = 1e8, 1.0
    param3d = Parameter.defaultUnit(1 / beta, rs, 3)
    param2d = Parameter.defaultUnit(1 / beta, rs, 2)

    @testset "Polarization: ZeroTemp vs. FiniteTemp" begin
        # for low temp, Π from zero temp and finite temp should be close
        qgrid = [-1.0, 0.0, 1e-16, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]
        ngrid = [n for n in -5:20]
        for (qi, q) in enumerate(qgrid)
            for (ni, n) in enumerate(ngrid)
                # println("q=$q, n=$n")
                @test isapprox(
                    Polarization.Polarization0_ZeroTemp(q, n, param3d),
                    Polarization.Polarization0_FiniteTemp(q, n, param3d),
                    rtol=1e-6
                )
            end
        end

        for (qi, q) in enumerate(qgrid)
            for (ni, n) in enumerate(ngrid)
                @test isapprox(
                    Polarization.Polarization0_ZeroTemp(q, n, param2d),
                    Polarization.Polarization0_FiniteTemp(q, n, param2d),
                    rtol=1e-4
                )
            end
        end
    end

    qgrid = [-1.0, 0.0, 1e-16, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]
    ngrid = [n for n in -5:20]

    @testset "Polarization0!: 3D ZeroTemp vs. FiniteTemp" begin
        # param = Parameter.defaultUnit(1 / beta, rs, 3)
        nmesh = MeshGrids.ImFreq(param3d.β, BOSON; grid=ngrid)
        PZ = MeshArray(nmesh, qgrid; dtype=Float64)
        PF = similar(PZ)

        Polarization.Polarization0_ZeroTemp!(PZ, param3d)
        Polarization.Polarization0_FiniteTemp!(PF, param3d)
        for ind in eachindex(PZ)
            @test isapprox(PZ[ind], PF[ind], rtol=1e-6)
        end
    end

    @testset "Polarization: wrapped" begin
        PZ = Polarization.Polarization0wrapped(100param3d.EF, 1e-10, qgrid, param3d, pifunc=Polarization.Polarization0_ZeroTemp!)
        PF = Polarization.Polarization0wrapped(100param3d.EF, 1e-10, qgrid, param3d, pifunc=Polarization.Polarization0_FiniteTemp!)

        PZ = Polarization.Polarization0wrapped(100param2d.EF, 1e-10, qgrid, param2d, pifunc=Polarization.Polarization0_ZeroTemp!)

        # PZ = Polarization.Polarization0wrapped(100param.EF, 1e-10, qgrid, param, pifunc=Polarization.Polarization0_ZeroTemp)
        # PF = Polarization.Polarization0wrapped(100param.EF, 1e-10, qgrid, param, pifunc=Polarization.Polarization0_FiniteTemp)
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
