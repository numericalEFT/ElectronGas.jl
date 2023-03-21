module Propogators

using ElectronGas
using ElectronGas.Interaction: RPAwrapped
using ElectronGas.Parameter

end

using Test

if abspath(PROGRAM_FILE) == @__FILE__
    using .Propogators

    @testset "Propogators" begin
        param = Propogators.Parameter.rydbergUnit(1, 1)
        rpai, rpad = Propogators.RPAwrapped(100 * param.EF, 1e-8, [1e-8, 0.5, 1.0, 2.0, 10.0], param)
    end

end