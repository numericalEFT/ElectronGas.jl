# provide and test containers of propagators for MC calculations.

module Propagators

using ElectronGas
using ElectronGas.Interaction: RPAwrapped
using ElectronGas.Parameter
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

export rpa, interaction, G0, initR, response

# Ri excludes the source term

# wrapper
struct Funcs{P,II,IT,RI,RT}
    param::P
    inti::II
    intt::IT
    Ri::RI
    Rt::RT
end

interaction(k, funcs::Funcs) = interaction(k, funcs.inti)
interaction(t, k, funcs::Funcs) = interaction(t, k, funcs.intt)
response(k, funcs::Funcs; norm=1) = response(k, funcs.Ri; norm=norm)
response(t, k, funcs::Funcs; norm=1) = response(t, k, funcs.Rt; norm=norm)
G0(t, k, funcs::Funcs) = G0(t, k, funcs.param)

# shift tau to [0, β)
function tau_fermi(t, β)
    if 0 <= t < β
        return t, 1.0
    elseif β <= t < 2β
        return t - β, -1.0
    elseif -β <= t < 0
        return t + β, -1.0
    else
        error("τ=$t out of range!")
    end
end

function tau_bose(t, β)
    if 0 <= t < β
        return t, 1.0
    elseif β <= t < 2β
        return t - β, 1.0
    elseif -β <= t < 0
        return t + β, 1.0
    else
        error("τ=$t out of range!")
    end
end

function rpa(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    Nlog = floor(Int, 2.0 * log10(param.β / mint))
    tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    rpad, rpai = Propagators.RPAwrapped(100 * param.EF, 1e-10, kgrid, param)
    rpadlr = Propagators.to_dlr(rpad)

    rpat = Propagators.dlr_to_imtime(rpadlr, tgrid)
    return rpai, rpat
end

function phonon(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    Nlog = floor(Int, 2.0 * log10(param.β / mint))
    tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    Euv = 100 * param.EF
    rtol = 1e-10
    β = param.β

    wn_mesh = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn = GreenFunc.MeshArray(1:2, wn_mesh, kgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(wn_mesh.grid)
            green_dyn[1, ni, ki], green_dyn[2, ni, ki] = Interaction.phonon(k, n, param)
        end
    end
    phonondlr = Propagators.to_dlr(green_dyn)
    phonont = Propagators.dlr_to_imtime(phonondlr, tgrid)

    return phonont
end

function interaction(k, prop; i=1)
    # return Interp.interp1D(view(prop.data, i, 1, :), prop.mesh[3], k)
    return Interp.linear1D(view(prop.data, i, 1, :), prop.mesh[3], k)
end

function interaction(t, k, prop; i=1)
    t, factor = tau_fermi(t, prop.mesh[2].β)
    # return factor * Interp.interpND(view(prop.data, i, :, :), prop.mesh[2:3], (t, k))
    return factor * Interp.linear2D(view(prop.data, i, :, :), prop.mesh[2], prop.mesh[3], t, k)
end

function G0(t, k, param)
    β = param.β
    t, factor = tau_fermi(t, β)
    ε = k^2 / 2 / param.me - param.μ
    f = 1 / (exp(-ε * β) + 1)
    result = factor * exp(-t * ε) * f
    if isnan(result) || isinf(result)
        println("t=$t, k=$k, ε=$ε, f=$f")
    end
    return result
end

function initR(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    Nlog = floor(Int, 2.0 * log10(param.β / mint))
    btgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)
    tgrid = GreenFunc.ImTime(param.β, true; symmetry=:pha, grid=btgrid)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    ri = GreenFunc.MeshArray(kgrid; dtype=ComplexF64)
    rt = GreenFunc.MeshArray(tgrid, kgrid; dtype=ComplexF64)
    ri[:] .= 0.0
    rt[:] .= 0.0

    return ri, rt
end

function R0(ri, rt, param)
    rdlr = Propagators.to_dlr(rt)
    rw = Propagators.to_imfreq(rdlr)
    kF = param.kF
    kgrid = rw.mesh[2]
    ikF = searchsortedfirst(kgrid, kF)
    return 1.0 + ri[ikF] + rw[1, ikF]
end

function response(k, ri; norm=1)
    # return Interp.interp1D(view(ri.data, :), ri.mesh[1], k)
    return 1.0 + Interp.linear1D(view(ri.data, :), ri.mesh[1], k) / norm
end

function response(t, k, rt; norm=1)
    t, factor = tau_fermi(t, rt.mesh[1].β)
    # return factor * Interp.interpND(view(rt.data, :, :), rt.mesh[:], (t, k))
    return factor * Interp.linear2D(view(rt.data, :, :), rt.mesh[1], rt.mesh[2], t, k) / norm
end

end

using Test, BenchmarkTools

if abspath(PROGRAM_FILE) == @__FILE__
    using .Propagators

    @testset "Propagators" begin

        param = Propagators.Parameter.defaultUnit(0.01, 1.0)

        rpai, rpat = Propagators.rpa(param)

        println(size(rpat))
        println(size(rpai))

        p = (0.00372, 0.0733)
        println(Propagators.interaction(p..., rpat))
        println(Propagators.interaction(p[2], rpai))
        @time Propagators.interaction(p..., rpat)
        @time Propagators.interaction(p[2], rpai)
        @time Propagators.G0(p..., param)

        @time Propagators.interaction(p..., rpat)
        @time Propagators.interaction(p[2], rpai)
        @time Propagators.G0(p..., param)

        Ri, Rt = Propagators.initR(param)
        @time Propagators.response(p..., Rt)
        @time Propagators.response(p[2], Ri)

        @time Propagators.response(p..., Rt)
        @time Propagators.response(p[2], Ri)

        pht = Propagators.phonon(param)
        @time Propagators.interaction(p..., pht)
        @time Propagators.interaction(p..., pht)
    end

end