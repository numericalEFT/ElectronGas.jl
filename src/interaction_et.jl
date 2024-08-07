
"""
    function ETdyson(V::Float64, F::Float64, Π::Float64)

Return (V - F) Π V / (1 - (V - F)Π), which is the dynamic part of the electron test-charge interaction W_et.

#Arguments:
- Vinv: inverse bare interaction
- F: Landau parameter
- Π: polarization
"""
function ETdyson(Vinv::Float64, F::Float64, Π::Float64)
    K = 0
    if Vinv ≈ Inf
        K = 0
    else
        K = Π / (Vinv / (1 - F * Vinv) - (Π)) / Vinv
    end
    @assert !isnan(K) "nan at Vinv=$Vinv, F=$F, Π=$Π"
    return K
end

"""
    function ETdysonreg(V::Float64, F::Float64, Π::Float64)

Return V Π / (1 - (V - F)Π), which is the dynamic part of the electron test-charge interaction W_et divided by (V - F).

#Arguments:
- Vinv: inverse bare interaction
- F: Landau parameter
- Π: polarization
"""
function ETdysonreg(Vinv::Float64, F::Float64, Π::Float64)
    K = 0
    if Vinv ≈ Inf
        K = 0
    else
        K = Π / (Vinv - (1 - F * Vinv) * Π)
    end
    @assert !isnan(K) "nan at Vinv=$Vinv, F=$F, Π=$Π"
    return K
end

function ETcorrection(q::Float64, n::Int, param;
    pifunc=Polarization0_ZeroTemp, landaufunc, Vinv_Bare=coulombinv, regular=false, massratio=1.0, kwargs...)
    Fs::Float64, Fa::Float64 = landaufunc(q, n, param; massratio=massratio, kwargs...)
    ETs::Float64, ETa::Float64 = 0.0, 0.0
    Vinvs::Float64, Vinva::Float64 = Vinv_Bare(q, param)
    @unpack spin = param

    if abs(q) < EPS
        q = EPS
    end

    Πs::Float64 = spin * pifunc(q, n, param; kwargs...) * massratio
    Πa::Float64 = spin * pifunc(q, n, param; kwargs...) * massratio
    if regular
        ETs = ETdysonreg(Vinvs, Fs, Πs)
        ETa = ETdysonreg(Vinva, Fa, Πa)
    else
        ETs = ETdyson(Vinvs, Fs, Πs)
        ETa = ETdyson(Vinva, Fa, Πa)
    end

    return ETs, ETa
end

"""
    function ET(q, n, param; pifunc = Polarization0_ZeroTemp, landaufunc = landauParameterTakada, Vinv_Bare = coulombinv, regular = false, kwargs...)

Dynamic part of electron test-charge interaction W_et. Returns the spin symmetric part and asymmetric part separately.

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
    \\frac{v_q^{\\pm} Π_0} {1 - (v_q^{\\pm} - f_q^{\\pm}) Π_0}.
```
otherwise, return
```math
    \\frac{(v_q^{\\pm} - f_q^{\\pm}) Π_0 v_q^{\\pm}} {1 - (v_q^{\\pm} - f_q^{\\pm}) Π_0}.
```
"""
function ET(q, n, param; pifunc=Polarization0_ZeroTemp, landaufunc, Vinv_Bare=coulombinv, regular=false, kwargs...)
    return ETcorrection(q, n, param; pifunc=pifunc, landaufunc=landaufunc, Vinv_Bare=Vinv_Bare, regular=regular, kwargs...)
end

"""
    function ETwrapped(Euv, rtol, sgrid::SGT, param; int_type=:ko,
        pifunc=Polarization0_ZeroTemp, landaufunc, Vinv_Bare=coulombinv, kwargs...) where {SGT}

Return dynamic part and instant part of ET-interaction Green's function separately. Each part is a MeshArray with inner state 
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
function ETwrapped(Euv, rtol, sgrid::SGT, param; int_type=:ko,
    pifunc=Polarization0_ZeroTemp, landaufunc, Vinv_Bare=coulombinv, kwargs...) where {SGT}
    # TODO: innerstate should be in the outermost layer of the loop. Hence, the functions such as ET and Vinv_Bare should be fixed with inner state as argument.
    @unpack β = param
    wn_mesh = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn = GreenFunc.MeshArray(1:2, wn_mesh, sgrid; dtype=ComplexF64)
    green_ins = GreenFunc.MeshArray(1:2, [0,], sgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(wn_mesh.grid)
            green_dyn[1, ni, ki], green_dyn[2, ni, ki] = ET(k, n, param; pifunc=pifunc, landaufunc=landaufunc, Vinv_Bare=Vinv_Bare, kwargs...)
        end
        green_ins[1, 1, ki], green_ins[2, 1, ki] = Vinv_Bare(k, param)
    end

    return green_dyn, green_ins
end
