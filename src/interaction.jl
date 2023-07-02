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
    pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, regular=false, massratio=1.0, kwargs...)
    Fs::Float64, Fa::Float64 = landaufunc(q, n, param; massratio=massratio, kwargs...)
    Ks::Float64, Ka::Float64 = 0.0, 0.0
    Vinvs::Float64, Vinva::Float64 = Vinv_Bare(q, param)
    @unpack spin = param

    if abs(q) < EPS
        q = EPS
    end

    Πs::Float64 = spin * pifunc(q, n, param; kwargs...) * massratio
    Πa::Float64 = spin * pifunc(q, n, param; kwargs...) * massratio
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
function RPA(q, n, param; pifunc=Polarization0_ZeroTemp, Vinv_Bare=coulombinv, regular=false, kwargs...)
    return bubblecorrection(q, n, param; pifunc=pifunc, landaufunc=landauParameter0, Vinv_Bare=Vinv_Bare, regular=regular, kwargs...)
end

"""
    function RPAwrapped(Euv, rtol, sgrid::SGT, param;
        pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, kwargs...) where {SGT}

Return dynamic part and instant part of RPA-interaction Green's function separately. Each part is a MeshArray with inner state 
(1: spin symmetric part, 2: asymmetric part), and ImFreq and q-grid mesh.

#Arguments:
 - Euv: Euv of DLRGrid
 - rtol: rtol of DLRGrid
 - sgrid: momentum grid
 - param: other system parameters
 - pifunc: caller to the polarization function
 - landaufunc: caller to the Landau parameter (exchange-correlation kernel)
 - Vinv_Bare: caller to the bare Coulomb interaction
"""
function RPAwrapped(Euv, rtol, sgrid::SGT, param;
    pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, kwargs...) where {SGT}
    # TODO: innerstate should be in the outermost layer of the loop. Hence, the functions such as RPA and Vinv_Bare should be fixed with inner state as argument.  
    @unpack β = param

    wn_mesh = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn = GreenFunc.MeshArray(1:2, wn_mesh, sgrid; dtype=ComplexF64)
    green_ins = GreenFunc.MeshArray(1:2, [0,], sgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(wn_mesh.grid)
            green_dyn[1, ni, ki], green_dyn[2, ni, ki] = RPA(k, n, param; pifunc=pifunc, Vinv_Bare=Vinv_Bare, regular=true, kwargs...)
        end
        green_ins[1, 1, ki], green_ins[2, 1, ki] = Vinv_Bare(k, param)
    end

    return green_dyn, green_ins
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
    A1 = (2.0 - C1 - C2) * (kF * ϵ0 * π^2) / (2 * e0^2 * me)
    A2 = (C2 - C1) * (kF * ϵ0 * π^2) / (2 * e0^2 * me)
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

"""
    function KOwrapped(Euv, rtol, sgrid::SGT, param; int_type=:ko,
        pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, kwargs...) where {SGT}

Return dynamic part and instant part of KO-interaction Green's function separately. Each part is a MeshArray with inner state 
(1: spin symmetric part, 2: asymmetric part), and ImFreq and q-grid mesh.

#Arguments:
 - Euv: Euv of DLRGrid
 - rtol: rtol of DLRGrid
 - sgrid: momentum grid
 - param: other system parameters
 - pifunc: caller to the polarization function
 - landaufunc: caller to the Landau parameter (exchange-correlation kernel)
 - Vinv_Bare: caller to the bare Coulomb interaction
"""
function KOwrapped(Euv, rtol, sgrid::SGT, param; int_type=:ko,
    pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, kwargs...) where {SGT}
    # TODO: innerstate should be in the outermost layer of the loop. Hence, the functions such as KO and Vinv_Bare should be fixed with inner state as argument.
    @unpack β = param
    wn_mesh = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn = GreenFunc.MeshArray(1:2, wn_mesh, sgrid; dtype=ComplexF64)
    green_ins = GreenFunc.MeshArray(1:2, [0,], sgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(wn_mesh.grid)
            green_dyn[1, ni, ki], green_dyn[2, ni, ki] = KO(k, n, param; pifunc=pifunc, landaufunc=landaufunc, Vinv_Bare=Vinv_Bare, kwargs...)
        end
        green_ins[1, 1, ki], green_ins[2, 1, ki] = Vinv_Bare(k, param)
    end

    return green_dyn, green_ins
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

"""
    function oldTmatrix(q, n, param; pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv, regular = false, kwargs...)

Dynamic part of the T-matrix. Returns the spin symmetric part and asymmetric part separately.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn (either one frequency or array of frequencies)
 - param: other system parameters
 - pifunc: caller to the polarization function 
 - landaufunc: caller to the Landau parameter (exchange-correlation kernel)
 - Vinv_Bare: caller to the bare Coulomb interaction
 - regular: regularized RPA or not

# Return:
```math
    \\frac{1} {\\frac{m}{4\\pi a_s} +\\Gamma(q, n)} - \\frac{4\\pi a_s}{m}.
```
"""
function oldTmatrix(q, n, param; ladderfunc=Polarization.Ladder0_FiniteTemp, kwargs...)
    as = 4π * param.as / param.me
    if param.as < 1e-6
        return as ./ (1 .+ as * ladderfunc(q, n, param; kwargs...))
    else
        return 1 ./ (1 / as .+ ladderfunc(q, n, param; kwargs...))
    end
end

"""
    function oldTmatrix_wrapped(Euv, rtol, sgrid::SGT, param; int_type=:ko,
        pifunc=Polarization0_ZeroTemp, landaufunc=landauParameterTakada, Vinv_Bare=coulombinv, kwargs...) where {SGT}

Return dynamic part and instant part of KO-interaction Green's function separately. Each part is a MeshArray with inner state 
(1: spin symmetric part, 2: asymmetric part), and ImFreq and q-grid mesh.

#Arguments:
 - Euv: Euv of DLRGrid
 - rtol: rtol of DLRGrid
 - sgrid: momentum grid
 - param: other system parameters
 - pifunc: caller to the polarization function
 - landaufunc: caller to the Landau parameter (exchange-correlation kernel)
 - Vinv_Bare: caller to the bare Coulomb interaction
"""
function oldTmatrix_wrapped(Euv, rtol, sgrid::SGT, param;
    ladderfunc=Polarization.Ladder0_FiniteTemp, kwargs...) where {SGT}
    # TODO: innerstate should be in the outermost layer of the loop. Hence, the functions such as KO and Vinv_Bare should be fixed with inner state as argument.
    @unpack β = param
    as = 4π * param.as / param.me
    wn_mesh_ph = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn_r = GreenFunc.MeshArray(wn_mesh_ph, sgrid; dtype=ComplexF64)
    green_tail_r = GreenFunc.MeshArray(wn_mesh_ph, sgrid; dtype=ComplexF64)

    wn_mesh_pha = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:pha)
    green_dyn_i = GreenFunc.MeshArray(wn_mesh_pha, sgrid; dtype=ComplexF64)
    green_tail_i = GreenFunc.MeshArray(wn_mesh_pha, sgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(wn_mesh_ph.grid)
            green_tail_r[ni, ki] = real.(Tmatrix(k, n, param; ladderfunc=Polarization.Ladder0_FreeElectron, kwargs...))
            green_dyn_r[ni, ki] = real.(Tmatrix(k, n, param; ladderfunc=ladderfunc, kwargs...)) - green_tail_r[ni, ki]
        end
        for (ni, n) in enumerate(wn_mesh_pha.grid)
            green_tail_i[ni, ki] = 1im .* imag.(Tmatrix(k, n, param; ladderfunc=Polarization.Ladder0_FreeElectron, kwargs...))
            green_dyn_i[ni, ki] = 1im .* imag.(Tmatrix(k, n, param; ladderfunc=ladderfunc, kwargs...)) - green_tail_i[ni, ki]
        end
    end

    return green_dyn_r, green_dyn_i, green_tail_r, green_tail_i
end


"""
    function T0(q::Float64, n::Int, param)

Vacuum two-particle T-matrix in matsubara frequency. 
#Arguments:
 - q: total incoming momentum
 - n: matsubara frequency given in integer s.t. Ωn=2πTn (either one frequency or array of frequencies)
 - param: other system parameters
"""

function T0matrix(q::Float64, n::Int, param)
    β, me = param.β, param.me
    as = 4π * param.as / me
    B = me^(3/2)/(4π)
    if as < 1e-6
        return as/(1. - as * B * sqrt(Complex(q^2/ (4 * me) - 2π * n/β * im)))
    else
        return 1/((1/as) - B * sqrt(Complex(q^2/ (4 * me) - 2π * n/β * im)))
    end
end

"""
    function Tmatrix(q::Float64, n::Int, param; ladderfunc=Polarization.Ladder0_FiniteTemp, vacuumladderfunc=Polarization.LadderVacuum0_FiniteTemp, kwargs...)

Many-body T-matrix in matsubara frequency. 
#Arguments:
 - q: total incoming momentum
 - n: matsubara frequency given in integer s.t. Ωn=2πTn (either one frequency or array of frequencies)
 - param: other system parameters
 - ladderfunc: Particle-particle bubble (ladder-diagram of G0G0)
 - vacuumladderfunc: Vacuum particle-particle bubble
"""

function Tmatrix(q::Float64, n::Int, param; ladderfunc=Polarization.Ladder0_FiniteTemp, vacuumladderfunc=Polarization.Ladder0Vacuum_FiniteTemp, kwargs...)
    as = 4π * param.as / param.me
    if param.as < 1e-6
        return as ./ (1 .+ as * (ladderfunc(q, n, param; kwargs...) - vacuumladderfunc(q, n, param; kwargs...)))
    else
        return 1 ./ (1 / as .+ (ladderfunc(q, n, param; kwargs...) - vacuumladderfunc(q, n, param; kwargs...)))
    end
end

"""
    function δTmatrix_wrapped(Euv, rtol, sgrid::SGT, param; ladderfunc=Polarization.Ladder0_FiniteTemp, kwargs...) where {SGT}

Difference between Many-body T-matrix and vacuum two-particle T-matrix in imaginary time. 
#Arguments:
- Euv: Euv of DLRGrid
- rtol: rtol of DLRGrid
- sgrid: momentum grid
- param: other system parameters
- ladderfunc: Particle-Particle bubble (ladder-diagram of G0G0)
"""

function δTmatrix_wrapped(Euv, rtol, sgrid::SGT, param; ladderfunc=Polarization.Ladder0_FiniteTemp, vacuumladderfunc=Polarization.Ladder0Vacuum_FiniteTemp, kwargs...) where {SGT}
    β= param.β
    Ωn_mesh_ph = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    δTmatrix_n = GreenFunc.MeshArray(Ωn_mesh_ph, sgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(Ωn_mesh_ph.grid)
            δTmatrix_n[ni, ki] = (Tmatrix(k, n, param; ladderfunc = ladderfunc, vacuumladderfunc = vacuumladderfunc) - T0matrix(k, n, param))
        end
    end
    δTmatrix_tau = δTmatrix_n |> to_dlr |> to_imtime
    return δTmatrix_tau
end

"""
    function T0maxtrix_imtime(tau::Float64, q::Float64, param)

Difference between Many-body T-matrix and vacuum two-particle T-matrix in imaginary time. 
#Arguments:
- tau: imaginary time
- q: total incoming momentum
- param: other system parameters
"""

function T0maxtrix_imtime(tau::Float64, q::Float64,  param)
    """ Two particle scattering in a vacuum T-matrix in imaginary time """
    β, me, as = param.β, param.me, param.as
    B = me^(3/2)/(4π)
    m = me * as^2

    # Pole Term contribution
    poletermnumerator = exp(tau/ m)
    poletermdenominator = 1- exp(β * (1. / (me*as^2) - q^2/ (4 * me)))
    poleterm = - poletermnumerator/ poletermdenominator
    # Integral Terms

    function integrand1(vars, q, tau, param)
        me, as = param.me, param.as
        m = me * as^2
        u = vars[1][1]
        y = u/(1-u)
        fraction = y^2/(y^2 + tau/ m)
        return 2 * exp(-y^2) * fraction
    end
    
    function integrand2(vars, q, tau, param)
        β, me = param.β, param.me
        u = vars[1][1]
        y = u/(1-u)
        factor = 1/(exp(β * (y^2/ tau + q^2/ (4 * me)))-1)
        return factor * integrand1(vars, q, tau, param)
    end

    tgrid = CompositeGrid.LogDensedGrid(
        :gauss,# The top layer grid is :gauss, optimized for integration. For interpolation use :cheb
        [0.0, 1.0],# The grid is defined on [0.0, 1.0] s.t. the integral in the original variable is from [0.0 , infinity]
        [0.0],# and is densed at 0.0 and β, as given by 2nd and 3rd parameter.
        5,# N of log grid
        0.005, # minimum interval length of log grid
        5 # N of bottom layer
    )
    data = [(integrand1(t, q, tau, param)+integrand2(t, q, tau, param)) for (ti, t) in  enumerate(tgrid.grid)]
    int_result = Interp.integrate1D(data, tgrid)
    integralfactor = 1/(B * pi * sqrt(tau))

    # Propagator factor
    propfactor = exp(-q^2/ (4 * me) * tau)

    return propfactor*(poleterm + int_result * integralfactor)
end

end
