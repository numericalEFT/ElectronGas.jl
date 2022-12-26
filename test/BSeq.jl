@testset "Bethe-Slapter equation" begin
    @testset "Cooper-pair linear response in 3D" begin
        dim, θ, rs = 3, 1e-2, 2.0
        param = Parameter.rydbergUnit(θ, rs, dim)
        channel = 0
        lamu, R_freq, F_freq = BSeq.linearResponse(param, channel)
        @test isapprox(lamu, -2.34540, rtol=1e-4)
        # lamu, R_freq = BSeq.linearResponse(param, 1)
        # @test isapprox(lamu, -1.23815, rtol=1e-5)

        # test stop condition. When SC, converge to 0.0
        dim, θ, rs = 3, 1 / 3200, 7.0
        param = Parameter.rydbergUnit(θ, rs, dim)
        channel = 0
        lamu, R_freq, F_freq = BSeq.linearResponse(param, channel)
        @test isapprox(lamu, 0.0, rtol=1e-10, atol=1e-10)

    end
    # @testset "Cooper-pair linear response in 2D" begin
    #     dim, θ, rs = 2, 1e-2, 1.5
    #     param = Parameter.rydbergUnit(θ, rs, dim)
    #     channel = 0
    #     lamu, R_freq = BSeq.linearResponse(param, channel)
    #     @test isapprox(lamu, -1.60555, rtol=1e-5)

    #     lamu, R_freq = BSeq.linearResponse(param, 1)
    #     @test isapprox(lamu, -1.05989, rtol=1e-5)
    # end
end