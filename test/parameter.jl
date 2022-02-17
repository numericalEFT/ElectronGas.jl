@testset "Parameter" begin
    # test constructor and reconstruct of Para
    
    beta = 1e4
    rs = 2.0
    param = Parameter.defaultUnit(1 / beta, rs)

    @test param.me == 0.5
    @test param.EF == 1.0
    @test param.kF == 1.0
    @test param.β == beta
    @test param.e0a == 0.0

    newbeta = 1e3
    espin = 1.0
    newparam = Parameter.DerivePara(param, β = newbeta, espin = espin)

    @test newparam.me == 0.5
    @test newparam.EF == 1.0
    @test newparam.kF == 1.0
    @test newparam.β == newbeta
    @test newparam.e0a == espin

end
