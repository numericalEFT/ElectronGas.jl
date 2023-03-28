# provide and test containers of propagators for MC calculations.

module Propagators

using ElectronGas
using ElectronGas.Interaction: RPAwrapped
using ElectronGas.Parameter
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

function rpa(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    Nlog = floor(Int, 2.0 * log10(param.beta / mint))
    tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)


    rpad, rpai = Propagators.RPAwrapped(100 * param.EF, 1e-10, kgrid, param)
    rpadlr = Propagators.to_dlr(rpad)

    rpat = Propagators.dlr_to_imtime(rpadlr, tgrid)
    return rpai, rpat
end

function interaction(k, prop; i=1)
    return Interp.interp1D(view(prop.data, i, 1, :), prop.mesh[3], k)
    # return Interp.linear1D(view(prop.data, i, 1, :), prop.mesh[3], k)
end

function interaction(t, k, prop; i=1)
    return Interp.interpND(view(prop.data, i, :, :), prop.mesh[2:3], (t, k))
    # return Interp.linear2D(view(prop.data, i, :, :), prop.mesh[2], prop.mesh[3], t, k)
end

end

using Test, BenchmarkTools

if abspath(PROGRAM_FILE) == @__FILE__
    using .Propagators

    @testset "Propagators" begin
        param = Propagators.Parameter.defaultUnit(0.01, 1.0)

        # mint = 0.001
        # Nlog = floor(Int, 2.0 * log10(param.beta / mint))
        # order = 6
        # tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

        # maxK = 10.0 * param.kF
        # minK = 0.001 * sqrt(param.T * param.me)
        # Nk = floor(Int, 2.0 * log10(maxK / minK))
        # kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)


        # rpad, rpai = Propagators.RPAwrapped(100 * param.EF, 1e-10, kgrid, param)
        # rpadlr = Propagators.to_dlr(rpad)

        # rpat = Propagators.dlr_to_imtime(rpadlr, tgrid)

        rpai, rpat = Propagators.rpa(param)

        println(size(rpat))
        println(size(rpai))

        p = (0.00372, 0.0733)
        println(Propagators.interaction(p..., rpat))
        println(Propagators.interaction(p[2], rpai))
        @time Propagators.interaction(p..., rpat)
        @time Propagators.interaction(p[2], rpai)

        @time Propagators.interaction(p..., rpat)
        @time Propagators.interaction(p[2], rpai)
    end

end