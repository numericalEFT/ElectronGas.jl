using ElectronGas
using CompositeGrids
using Test

@testset "ElectronGas.jl" begin
    # Write your tests here.
    if isempty(ARGS)
        include("parameter.jl")
        include("polarization.jl")
        include("interaction.jl")
        include("legendreinteraction.jl")
        # include("selfenergy.jl")
    else
        include(ARGS[1])
    end

end
