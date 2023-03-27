# provide and test containers of propagators for MC calculations.

module Propagators

using ElectronGas
using ElectronGas.Interaction: RPAwrapped
using ElectronGas.Parameter
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

end

using Test

if abspath(PROGRAM_FILE) == @__FILE__
    using .Propagators

    @testset "Propagators" begin
        param = Propagators.Parameter.defaultUnit(0.01, 1.0)

        mint = 0.001
        Nlog = floor(Int, 2.0 * log10(param.beta / mint))
        order = 6
        tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

        maxK = 10.0 * param.kF
        minK = 0.001 * sqrt(param.T * param.me)
        Nk = floor(Int, 2.0 * log10(maxK / minK))
        kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)


        rpad, rpai = Propagators.RPAwrapped(100 * param.EF, 1e-10, kgrid, param)
        rpadlr = Propagators.to_dlr(rpad)

        rpat = Propagators.dlr_to_imtime(rpadlr, tgrid)
        println(size(rpat))

        p = (0.00372, 0.0733)
        println(Propagators.Interp.interpND(rpat.data[1, :, :], (tgrid, kgrid), p))
        println(Propagators.Interp.linearND(rpat.data[1, :, :], (tgrid, kgrid), p))
    end

end