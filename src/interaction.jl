module Interaction

# using Parameters, GreenFunc
# include(srcdir*"/parameter.jl")
# using .Parameter
# include(srcdir*"/convention.jl")
# using .Convention
# include(srcdir*"/polarization.jl")
# using .Polarization

using ..Parameter, ..Convention, ..Polarization
using ..Parameters, ..CompositeGrids, ..GreenFunc

export RPA, KO, RPAwrapped, KOwrapped, coulomb

# if !@isdefined Para
#     include(rundir*"/para.jl")
#     using .Para
# end

# println(Parameter.Param)
# @unpack me, kF, rs, e0, β , Λs, ϵ0= Parameter.Param

function inf_sum(q, n)
    # Calculate a series sum for Takada anzats
    # See Takada(doi:10.1103/PhysRevB.47.5202)(Eq.2.16).
    a = q * q
    sum = 1.0
    i = 0
    j = 1.0
    k = 2.0
    for i in 1:n
        sum = sum + a / j / k
        a = a * q * q
        j = j * (i + 1.0)
        k = k * (i + 2.0)
    end
    return 1.0 / sum / sum
end

"""
    function coulomb(q,param)

Bare interaction in momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function coulomb(q, param)
    @unpack me, kF, rs, e0s, e0a, β, Λs, Λa, ϵ0 = param
    if e0s ≈ 0.0
        Vs = 0.0
    else
        if (q^2 + Λs) ≈ 0.0
            Vs = Inf
        else
            Vs = e0s^2 / ϵ0 / (q^2 + Λs)
        end
    end
    if e0a ≈ 0.0
        Va = 0.0
    else
        if (q^2 + Λa) ≈ 0.0
            Va = Inf
        else
            Va = e0a^2 / ϵ0 / (q^2 + Λa)
        end
    end
    return Vs, Va
end

"""
    function coulombinv(q,param)

Inverse of bare interaction in momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function coulombinv(q, param)
    @unpack me, kF, rs, e0s, e0a, β, Λs, Λa, ϵ0 = param
    if e0s^2 ≈ 0.0
        Vinvs = Inf
    else
        Vinvs = ϵ0 * (q^2 + Λs) / e0s^2
    end
    if e0a^2 ≈ 0.0
        Vinva = Inf
    else
        Vinva = ϵ0 * (q^2 + Λa) / e0a^2
    end
    return  Vinvs, Vinva
end

"""
    function bubbledyson(Vinv::Float64, F::Float64, Π::Float64)

Return (V - F)^2 Π / (1 - (V - F)Π), which is the dynamic part of effective interaction.

#Arguments:
- Vinv: inverse bare interaction
- F: Landau parameter
- Π: polarization
"""
function bubbledyson(Vinv::Float64, F::Float64, Π::Float64)
    K = 0
    if Vinv ≈ Inf
        if F ≈ 0
            K = 0
        else
            K = Π / ( 1.0 / (-F) - (Π)) * (-F)
        end
    else
        K = Π / (Vinv / (1 - F * Vinv) - (Π)) * (1 - F * Vinv) / Vinv
    end
    @assert !isnan(K) "nan at Vinv=$Vinv, F=$F, Π=$Π"
    return K
end

"""
    function bubbledysonreg(Vinv::Float64, F::Float64, Π::Float64)

Return (V - F) Π / (1 - (V - F)Π), which is the dynamic part of effective interaction divided by (V - F).

#Arguments:
- Vinv: inverse bare interaction
- F: Landau parameter
- Π: polarization
"""
function bubbledysonreg(Vinv::Float64, F::Float64, Π::Float64)
    K = 0
    if Vinv ≈ Inf
        if F ≈ 0
            K = 0
        else
            K = Π / ( 1.0 / (-F) - (Π))
        end
    else
        K = Π / (Vinv / (1 - F * Vinv) - (Π))
    end
    @assert !isnan(K) "nan at Vinv=$Vinv, F=$F, Π=$Π"
    return K
end

function bubblecorrection(q::Float64, n::Int, param;
    pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv, isregularized = false)
    Fs::Float64, Fa::Float64 = landaufunc(q, n, param)
    Ks::Float64, Ka::Float64 = 0.0, 0.0
    # Vs::Float64, Va::Float64 = V_Bare(q, param)
    Vinvs::Float64, Vinva::Float64 = Vinv_Bare(q, param)
    @unpack spin = param

    if abs(q) > EPS
        Π::Float64 = spin * pifunc(q, n, param)
        if isregularized
            Ks = bubbledysonreg(Vinvs, Fs, Π)
            Ka = bubbledysonreg(Vinva, Fa, Π)
        else
            Ks = bubbledyson(Vinvs, Fs, Π)
            Ka = bubbledyson(Vinva, Fa, Π)
        end
    else
        Ks, Ka = 0.0, 0.0
    end

    return Ks, Ka
end

"""
    function RPA(q, n, param)

Dynamic part of RPA interaction, with polarization approximated by zero temperature Π0.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function RPA(q, n, param; pifunc = Polarization0_ZeroTemp, Vinv_Bare = coulombinv, isregularized = false)
    return bubblecorrection(q, n, param; pifunc = pifunc, landaufunc = landauParameter0, Vinv_Bare = Vinv_Bare, isregularized = isregularized)
end

function RPAwrapped(Euv, rtol, sgrid::SGT, param;
    pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv) where {SGT}

    @unpack β = param
    gs = GreenFunc.Green2DLR{Float64}(:rpa, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry = :ph, rtol = rtol)
    ga = GreenFunc.Green2DLR{Float64}(:rpa, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry = :ph, rtol = rtol)
    green_dyn_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size, gs.timeGrid.size))
    green_ins_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size))
    green_dyn_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size, ga.timeGrid.size))
    green_ins_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(gs.dlrGrid.n)
            green_dyn_s[1, 1, ki, ni], green_dyn_a[1, 1, ki, ni] = RPA(k, n, param; pifunc = pifunc, Vinv_Bare = Vinv_Bare)
        end
        green_ins_s[1, 1, ki], green_ins_a[1, 1, ki] = Vinv_Bare(k, param)
    end
    gs.dynamic = green_dyn_s
    gs.instant = green_ins_s
    ga.dynamic = green_dyn_a
    ga.instant = green_ins_a
    return gs, ga
end

"""
    function landauParameterTakada(q, n, param)

G factor with Takada's anzats. See Takada(doi:10.1103/PhysRevB.47.5202)(Eq.2.13-2.16).
Now Landau parameter F. F(+-)=G(+-)*V

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function landauParameterTakada(q, n, param)
    @unpack me, kF, rs, e0s, e0, e0a, β, Λs, Λa, ϵ0 = param
    if e0 ≈ 0.0
        return 0.0, 0.0
    end
    r_s_dl = sqrt(4 * 0.521 * rs / π)
    C1 = 1 - r_s_dl * r_s_dl / 4.0 * (1 + 0.07671 * r_s_dl * r_s_dl * ((1 + 12.05 * r_s_dl) * (1 + 12.05 * r_s_dl) + 4.0 * 4.254 / 3.0 * r_s_dl * r_s_dl * (1 + 7.0 / 8.0 * 12.05 * r_s_dl) + 1.5 * 1.363 * r_s_dl * r_s_dl * r_s_dl * (1 + 8.0 / 9.0 * 12.05 * r_s_dl)) / (1 + 12.05 * r_s_dl + 4.254 * r_s_dl * r_s_dl + 1.363 * r_s_dl * r_s_dl * r_s_dl) / (1 + 12.05 * r_s_dl + 4.254 * r_s_dl * r_s_dl + 1.363 * r_s_dl * r_s_dl * r_s_dl))
    C2 = 1 - r_s_dl * r_s_dl / 4.0 * (1 + r_s_dl * r_s_dl / 8.0 * (log(r_s_dl * r_s_dl / (r_s_dl * r_s_dl + 0.990)) - (1.122 + 1.222 * r_s_dl * r_s_dl) / (1 + 0.533 * r_s_dl * r_s_dl + 0.184 * r_s_dl * r_s_dl * r_s_dl * r_s_dl)))
    D = inf_sum(r_s_dl, 100)
    A1 = (2.0 - C1 - C2) / 4.0 / e0^2 * π
    A2 = (C2 - C1) / 4.0 / e0^2 * π
    B1 = 6 * A1 / (D + 1.0)
    B2 = 2 * A2 / (1.0 - D)
    F_s = A1 * e0^2 / ϵ0 / (1.0 + B1 * q^2) + A2 * e0^2 / ϵ0 / (1.0 + B2 * q^2)
    F_a = A1 * e0^2 / ϵ0 / (1.0 + B1 * q^2) - A2 * e0^2 / ϵ0 / (1.0 + B2 * q^2)
    return F_s, F_a
end

@inline function landauParameter0(q, n, param)
    return 0.0, 0.0
end

"""
    function KO(q, n, param)

Dynamic part of KO interaction, with polarization approximated by zero temperature Π0.
Returns the spin symmetric part and asymmetric part separately.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function KO(q, n, param; pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv, isregularized = false)
    return bubblecorrection(q, n, param; pifunc = pifunc, landaufunc = landaufunc, Vinv_Bare = coulombinv, isregularized = isregularized)
end

function KOwrapped(Euv, rtol, sgrid::SGT, param;
    pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv) where {SGT}

    @unpack β = param
    gs = GreenFunc.Green2DLR{Float64}(:ko, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry = :ph, rtol = rtol)
    ga = GreenFunc.Green2DLR{Float64}(:ko, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry = :ph, rtol = rtol)
    green_dyn_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size, gs.timeGrid.size))
    green_ins_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size))
    green_dyn_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size, ga.timeGrid.size))
    green_ins_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(gs.dlrGrid.n)
            green_dyn_s[1, 1, ki, ni], green_dyn_a[1, 1, ki, ni] = KO(k, n, param; pifunc = pifunc, landaufunc = landaufunc, Vinv_Bare = Vinv_Bare)
        end
        green_ins_s[1, 1, ki], green_ins_a[1, 1, ki] = Vinv_Bare(k, param)
    end
    gs.dynamic = green_dyn_s
    gs.instant = green_ins_s
    ga.dynamic = green_dyn_a
    ga.instant = green_ins_a
    return gs, ga
end

end
