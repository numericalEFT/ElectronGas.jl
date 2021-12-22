using ElectronGas
using CompositeGrids
using Test

@testset "ElectronGas.jl" begin
    # Write your tests here.
    if isempty(ARGS)
        include("polarization.jl")
    else
        include(ARGS[1])
    end

end
