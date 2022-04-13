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

export RPA, KO, RPAwrapped, KOwrapped, coulomb, coulomb_2d, landauParameterMoroni

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
    @unpack me, kF, rs, e0s, e0a, β, Λs, Λa, ϵ0, gs, ga = param
    if gs ≈ 0.0
        Vs = 0.0
    else
        if (q^2 + Λs) ≈ 0.0
            Vs = Inf
        else
            Vs = e0s^2 / ϵ0 / (q^2 + Λs) * gs
        end
    end
    if ga ≈ 0.0
        Va = 0.0
    else
        if (q^2 + Λa) ≈ 0.0
            Va = Inf
        else
            Va = e0a^2 / ϵ0 / (q^2 + Λa) * ga
        end
    end
    return Vs, Va
end

"""
    function coulomb_2d(q,param)

Bare interaction in 2D momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function coulomb_2d(q, param)
    @unpack me, kF, rs, e0s, e0a, β, Λs, Λa, ϵ0, gs, ga = param
    if gs ≈ 0.0
        Vs = 0.0
    else
        if (q^2 + Λs) ≈ 0.0
            Vs = Inf
        else
            Vs = e0s^2 / 2ϵ0 / √(q^2 + Λs) * gs
        end
    end
    if ga ≈ 0.0
        Va = 0.0
    else
        if (q^2 + Λa) ≈ 0.0
            Va = Inf
        else
            # Va = e0a^2 / 2ϵ0 / √(q^2 + Λa)
            Va = e0a^2 / ϵ0 / (q^2 + Λa) * ga
        end
    end
    return Vs, Va
end

"""
    function coulombinv(q,param)

Inverse of bare interaction in 3D momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function coulombinv(q, param)
    @unpack me, kF, rs, e0s, e0a, β, Λs, Λa, ϵ0, gs, ga = param
    if gs ≈ 0.0
        Vinvs = Inf
    else
        Vinvs = ϵ0 * (q^2 + Λs) / e0s^2 / gs
    end
    if ga ≈ 0.0
        Vinva = Inf
    else
        Vinva = ϵ0 * (q^2 + Λa) / e0a^2 / ga
    end
    return Vinvs, Vinva
end

"""
    function coulombinv_2d(q,param)

Inverse of bare interaction in 2D momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function coulombinv_2d(q, param)
    @unpack me, kF, rs, e0s, e0a, β, Λs, Λa, ϵ0, gs, ga = param
    if gs ≈ 0.0
        Vinvs = Inf
    else
        Vinvs = 2ϵ0 * √(q^2 + Λs) / e0s^2 / gs
    end
    if ga ≈ 0.0
        Vinva = Inf
    else
        Vinva = ϵ0 * (q^2 + Λa) / e0a^2 / ga
        # Vinva = 2ϵ0 * √(q^2 + Λa) / e0a^2
    end
    return Vinvs, Vinva
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
            K = Π / (1.0 / (-F) - (Π)) * (-F)
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
            K = Π / (1.0 / (-F) - (Π))
        end
    else
        K = Π / (Vinv / (1 - F * Vinv) - (Π))
    end
    @assert !isnan(K) "nan at Vinv=$Vinv, F=$F, Π=$Π"
    return K
end

function bubblecorrection(q::Float64, n::Int, param;
    pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, regular=false, kwargs...)
    Fs::Float64, Fa::Float64 = landaufunc(q, n, param; kwargs...)
    Ks::Float64, Ka::Float64 = 0.0, 0.0
    Vinvs::Float64, Vinva::Float64 = Vinv_Bare(q, param)
    @unpack spin = param

    if abs(q) < EPS
        q = EPS
    end

    Πs::Float64 = spin * pifunc(q, n, param)
    Πa::Float64 = spin * pifunc(q, n, param)
    if regular
        Ks = bubbledysonreg(Vinvs, Fs, Πs)
        Ka = bubbledysonreg(Vinva, Fa, Πa)
    else
        Ks = bubbledyson(Vinvs, Fs, Πs)
        Ka = bubbledyson(Vinva, Fa, Πa)
    end

    return Ks, Ka
end

"""
    function RPA(q, n, param; pifunc = Polarization0_ZeroTemp, Vinv_Bare = coulombinv, regular = false)

    Dynamic part of RPA interaction. 

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
 - pifunc: caller to the polarization function 
 - Vinv_Bare: caller to the bare Coulomb interaction
 - regular: regularized RPA or not

# Return:
If set to be regularized, it returns the dynamic part of effective interaction divided by ``v_q - f_q``
```math
    \\frac{v_q^{\\pm} Π_0} {1 - v_q^{\\pm} Π_0}.
```
otherwise, return
```math
    \\frac{(v_q^{\\pm})^2 Π_0} {1 - v_q^{\\pm} Π_0}.
```
"""
function RPA(q, n, param; pifunc=Polarization0_ZeroTemp, Vinv_Bare=coulombinv, regular=false)
    return bubblecorrection(q, n, param; pifunc=pifunc, landaufunc=landauParameter0, Vinv_Bare=Vinv_Bare, regular=regular)
end

function RPAwrapped(Euv, rtol, sgrid::SGT, param;
    pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv) where {SGT}

    @unpack β = param
    green_s = GreenFunc.Green2DLR{Float64}(:rpa, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry=:ph, rtol=rtol)
    green_a = GreenFunc.Green2DLR{Float64}(:rpa, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry=:ph, rtol=rtol)
    green_dyn_s = zeros(Float64, (green_s.color, green_s.color, green_s.spaceGrid.size, green_s.timeGrid.size))
    green_ins_s = zeros(Float64, (green_s.color, green_s.color, green_s.spaceGrid.size))
    green_dyn_a = zeros(Float64, (green_a.color, green_a.color, green_a.spaceGrid.size, green_a.timeGrid.size))
    green_ins_a = zeros(Float64, (green_a.color, green_a.color, green_a.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green_s.dlrGrid.n)
            green_dyn_s[1, 1, ki, ni], green_dyn_a[1, 1, ki, ni] = RPA(k, n, param; pifunc=pifunc, Vinv_Bare=Vinv_Bare, regular=true)
        end
        green_ins_s[1, 1, ki], green_ins_a[1, 1, ki] = Vinv_Bare(k, param)
    end
    green_s.dynamic = green_dyn_s
    green_s.instant = green_ins_s
    green_a.dynamic = green_dyn_a
    green_a.instant = green_ins_a
    return green_s, green_a
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
function landauParameterTakada(q, n, param; kwargs...)
    @unpack me, kF, rs, e0s, e0, e0a, β, Λs, Λa, ϵ0 = param
    if e0 ≈ 0.0
        return 0.0, 0.0
    end
    r_s_dl = sqrt(4 * 0.521 * rs / π)
    C1 = 1 - r_s_dl * r_s_dl / 4.0 * (1 + 0.07671 * r_s_dl * r_s_dl * ((1 + 12.05 * r_s_dl) * (1 + 12.05 * r_s_dl) + 4.0 * 4.254 / 3.0 * r_s_dl * r_s_dl * (1 + 7.0 / 8.0 * 12.05 * r_s_dl) + 1.5 * 1.363 * r_s_dl * r_s_dl * r_s_dl * (1 + 8.0 / 9.0 * 12.05 * r_s_dl)) / (1 + 12.05 * r_s_dl + 4.254 * r_s_dl * r_s_dl + 1.363 * r_s_dl * r_s_dl * r_s_dl) / (1 + 12.05 * r_s_dl + 4.254 * r_s_dl * r_s_dl + 1.363 * r_s_dl * r_s_dl * r_s_dl))
    C2 = 1 - r_s_dl * r_s_dl / 4.0 * (1 + r_s_dl * r_s_dl / 8.0 * (log(r_s_dl * r_s_dl / (r_s_dl * r_s_dl + 0.990)) - (1.122 + 1.222 * r_s_dl * r_s_dl) / (1 + 0.533 * r_s_dl * r_s_dl + 0.184 * r_s_dl * r_s_dl * r_s_dl * r_s_dl)))
    D = inf_sum(r_s_dl, 100)
    #A1 = (2.0 - C1 - C2) / 4.0 / e0^2 * π
    A1 = (2.0-C1-C2)*(kF*ϵ0*π^2)/(2*e0^2*me)
    A2 = (C2 - C1)*(kF*ϵ0*π^2)/(2*e0^2*me)
    B1 = 6 * A1 / (D + 1.0)
    B2 = 2 * A2 / (1.0 - D)
    F_s = A1 * e0^2 / ϵ0 / (kF^2 + B1 * q^2) + A2 * e0^2 / ϵ0 / (kF^2 + B2 * q^2)
    F_a = A1 * e0^2 / ϵ0 / (kF^2 + B1 * q^2) - A2 * e0^2 / ϵ0 / (kF^2 + B2 * q^2)
    return F_s, F_a
end

"""
    function landauParameterMoroni(q, n, param)

Analytic expression of G+ factor from diffusion Monte Carlo.(doi:10.1103/PhysRevB.57.14569).
Return Landau parameter F+=(G+)*V

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function landauParameterMoroni(q, n, param; kwargs...)
    @unpack me, kF, rs, e0s, e0, e0a, β, ϵ0, a0 = param

    n0 = (kF)^3 / 3 / π^2
    if e0 ≈ 0.0
        return 0.0, 0.0
    end
    xc = -0.10498
    bc = 3.72744
    cc = 12.9352
    Ac = 0.0621814
    function E_corr(rs_dum)
        result = 0
        x = sqrt(rs_dum)
        X = (x^2 + bc * x + cc)
        X0 = (xc^2 + bc * xc + cc)
        Q = sqrt(4 * cc - bc^2)
        #print("x^2/X=",x^2/X,"\n")
        #print("x-xc^2/X=",(x-xc)^2/X, "\n")
        result = Ac * (log(x^2 / X) + 2 * bc / Q * atan(Q / (2 * x + bc)) - bc * xc / X0 * (log((x - xc)^2 / X) + 2 * (bc + 2 * xc) / Q * atan(Q / (2 * x + bc))))
        # This expression is in Rydberg unit
        result = result / (2 * me * a0^2)
        return result
    end

    # step for numerical derivatives 
    step = 0.0000001
    x = sqrt(rs)
    # Calculate parameter A

    rs1 = (n0 / (n0 + step))^(1.0 / 3.0) * rs
    rs_1 = (n0 / (n0 - step))^(1.0 / 3.0) * rs
    deriv_1 = (E_corr(rs1) - E_corr(rs)) / step
    deriv_2 = (E_corr(rs1) + E_corr(rs_1) - 2 * E_corr(rs)) / step^2
    A = 0.25 - kF^2 / 4 / π / e0^2 * (2 * deriv_1 + n0 * deriv_2)
    #println("A=$(A)")

    # Calculate parameter B
    a1 = 2.15
    a2 = 0.435
    b1 = 1.57
    b2 = 0.409
    B = (1 + a1 * x + a2 * x^3) / (3 + b1 * x + b2 * x^3)
    #println("B=$(B)")

    # Calculate parameter C
    step = 0.0000001
    deriv_1 = (E_corr(rs + step) - E_corr(rs)) / step
    #deriv_1 = E_corr_p(rs)
    C = -π / 2 / kF / e0^2 * (E_corr(rs) + rs * deriv_1)
    #println("C=$(C)")

    D = B / (A - C)
    α = 1.5 / rs^0.25 * A / B / D
    β_0 = 1.2 / B / D
    Q = q / kF
    G_s = C * Q^2 + B * Q^2 / (D + Q^2) + α * Q^4 * exp(-β_0 * Q^2)
    F_s = 4 * π * e0^2 * G_s / q^2
    return F_s
end


@inline function landauParameter0(q, n, param; kwargs...)
    return 0.0, 0.0
end

@inline function landauParameterConst(q, n, param; Fs=0.0, Fa=0.0, massratio=1.0, kwargs...)
    return Fs / param.NF / massratio, Fa / param.NF / massratio
end

@inline function counterterm(q, n, param; landaufunc, kwargs...)
    fs, fa = landaufunc(q, n, param; kwargs...)
    return fs, fa
end

@inline function countertermConst(q, n, param; landaufunc, Cs=0.0, Ca=0.0, massratio=1.0, kwargs...)
    return Cs / param.NF / massratio, Ca / param.NF / massratio
end


"""
    function KO(q, n, param; pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv, regular = false, kwargs...)

Dynamic part of KO interaction. Returns the spin symmetric part and asymmetric part separately.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
 - pifunc: caller to the polarization function 
 - landaufunc: caller to the Landau parameter (exchange-correlation kernel)
 - Vinv_Bare: caller to the bare Coulomb interaction
 - regular: regularized RPA or not

# Return:
If set to be regularized, it returns the dynamic part of effective interaction divided by ``v_q - f_q``
```math
    \\frac{(v_q^{\\pm} - f_q^{\\pm}) Π_0} {1 - (v_q^{\\pm} - f_q^{\\pm}) Π_0}.
```
otherwise, return
```math
    \\frac{(v_q^{\\pm} - f_q^{\\pm})^2 Π_0} {1 - (v_q^{\\pm} - f_q^{\\pm}) Π_0}.
```
"""
function KO(q, n, param; pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, regular=false, kwargs...)
    return bubblecorrection(q, n, param; pifunc=pifunc, landaufunc=landaufunc, Vinv_Bare=Vinv_Bare, regular=regular, kwargs...)
end

function KOwrapped(Euv, rtol, sgrid::SGT, param; int_type=:ko,
    pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, kwargs...) where {SGT}

    @unpack β = param
    green_s = GreenFunc.Green2DLR{Float64}(int_type, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry=:ph, rtol=rtol)
    green_a = GreenFunc.Green2DLR{Float64}(int_type, GreenFunc.IMFREQ, β, false, Euv, sgrid, 1; timeSymmetry=:ph, rtol=rtol)
    green_dyn_s = zeros(Float64, (green_s.color, green_s.color, green_s.spaceGrid.size, green_s.timeGrid.size))
    green_ins_s = zeros(Float64, (green_s.color, green_s.color, green_s.spaceGrid.size))
    green_dyn_a = zeros(Float64, (green_a.color, green_a.color, green_a.spaceGrid.size, green_a.timeGrid.size))
    green_ins_a = zeros(Float64, (green_a.color, green_a.color, green_a.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green_s.dlrGrid.n)
            green_dyn_s[1, 1, ki, ni], green_dyn_a[1, 1, ki, ni] = KO(k, n, param; pifunc=pifunc, landaufunc=landaufunc, Vinv_Bare=Vinv_Bare, kwargs...)
        end
        green_ins_s[1, 1, ki], green_ins_a[1, 1, ki] = Vinv_Bare(k, param)
    end
    green_s.dynamic = green_dyn_s
    green_s.instant = green_ins_s
    green_a.dynamic = green_dyn_a
    green_a.instant = green_ins_a
    return green_s, green_a
end

"""
    function KO_total(q, n, param; pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv, counter_term = counterterm, kwargs...)

Dynamic part of KO interaction. Returns the spin symmetric part and asymmetric part separately.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
 - pifunc: caller to the polarization function 
 - landaufunc: caller to the Landau parameter (exchange-correlation kernel)
 - Vinv_Bare: caller to the bare Coulomb interaction
 - counter_term: counterterm, by default, it is the landaufunc

# Return:
Return the total effective interaction
```math
    W^{\\pm} = \\frac{(v_q^{\\pm} - f_q^{\\pm}) Π_0} {1 - (v_q^{\\pm} - f_q^{\\pm}) Π_0} + C_q^{\\pm}.
```
which reduces to the convential KO interaction if ``C_q^{\\pm} \\equiv f_q^{\\pm}``
"""
function KO_total(q, n, param; pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, counter_term=counterterm, kwargs...)
    @unpack spin = param

    if abs(q) < EPS
        q = EPS
    end

    Π::Float64 = spin * pifunc(q, n, param; kwargs...)

    fs, fa = landaufunc(q, n, param; kwargs...)
    Cs, Ca = counter_term(q, n, param; landaufunc=landaufunc, kwargs...)
    Vinvs, Vinva = Vinv_Bare(q, param)

    if param.gs ≈ 0.0
        Ka = (-fs) / (1 - (-fs) * Π) + Cs
    else
        Ks = 1.0 / (Vinvs / (1 - fs * Vinvs) - Π) + Cs
    end
    if param.ga ≈ 0.0
        Ka = (-fa) / (1 - (-fa) * Π) + Ca
    else
        Ka = 1.0 / (Vinva / (1 - fa * Vinva) - Π) + Ca
    end
    return Ks, Ka
end

function RPA_total(q, n, param; pifunc=Polarization0_ZeroTemp, Vinv_Bare=coulombinv, kwargs...)
    @unpack spin = param

    if abs(q) < EPS
        q = EPS
    end

    Π::Float64 = spin * pifunc(q, n, param)
    Vinvs, Vinva = Vinv_Bare(q, param)

    Ws = 1.0 / (Vinvs - (Π))
    Wa = 1.0 / (Vinva - (Π))
    return Ws, Wa
end

end
