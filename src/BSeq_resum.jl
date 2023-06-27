"""
Bethe-Slapter-type equation solver and the application to Cooper-pair linear response approach. 
"""
module BSeq_resum

using ..Parameter, ..Convention, ..LegendreInteraction, ..Interaction
using ..Parameters, ..GreenFunc, ..Lehmann, ..CompositeGrids
using ..SelfEnergy
using ..JLD2
using ..BSeq

const freq_sep = 0.01

"""
    function initFR(Euv, rtol, sgrid, param)

Initalize the Bethe-Slapter amplitude `F` in imaginary-frequency space and `R ≡ F(GG)⁻¹` in imaginary-time space.

# Arguments:
- `Euv`: the UV energy scale of the spectral density. parameter for DLR grids.
- `rtol`: tolerance absolute error. parameter for DLR grids.
- `sgrid`: momentum grid of F and R
- `param`: parameters of ElectronGas.

# Return
- ``F(\\omega_n, k)=0`` as a `GreenFunc.MeshArray`
- the dynamical part of `R` in the imaginary-time space as a `GreenFunc.MeshArray`. In the imaginary-frequency space,
```math
    R(\\omega_n, k) =  \\frac{1}{\\Omega_c^2+e^2}\\left(1- \\frac{2\\omega_n^2}{\\omega_n^2+\\Omega_c^2} \\right)
```
where ``e=k^2/(2m)-\\mu`` and ``\\Omega_c =0.01``.
- the instant part of `R` as a `GreenFunc.MeshArray`, ``R_{\\mathrm{ins}}(k)=0``.
"""
# function initFR(Euv, rtol, sgrid, param)
#     @unpack β, me, μ = param
#     Ω_c = freq_sep

#     wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
#     R_freq = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
#     F_freq = similar(R_freq)
#     for ind in eachindex(R_freq)
#         ni, ki = ind[1], ind[2]
#         ωn, k = wn_mesh[ni], sgrid[ki]
#         e = k^2 / 2 / me - μ
#         R_freq[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
#     end
#     R_imt = real(R_freq |> to_dlr |> to_imtime)
#     R_ins = GreenFunc.MeshArray([1], sgrid; dtype=Float64, data=zeros(1, sgrid.size))
#     return F_freq, R_imt, R_ins
# end

function initFR_resum(Euv, rtol, sgrid, param)
    @unpack β, me, μ = param
    Ω_c = freq_sep

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    R_freq = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
    F_freq = similar(R_freq)
    for ind in eachindex(R_freq)
        ni, ki = ind[1], ind[2]
        ωn, k = wn_mesh[ni], sgrid[ki]
        e = k^2 / 2 / me - μ
        R_freq[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
    end
    R_imt = real(R_freq |> to_dlr |> to_imtime)
    R_ins = GreenFunc.MeshArray([1], sgrid; dtype=Float64, data=zeros(1, sgrid.size))
    F_fs = GreenFunc.MeshArray(wn_mesh, [1,]; dtype=Float64, data=zeros(length(wn_mesh), 1))
    return F_freq, F_fs, R_imt, R_ins
end

"""
    function initFR(Euv, rtol, param)

Initalize the Uniform amplitude `R` and zero amplitude 'F' in imaginary-frequency space. Both function have no momentum dependence. 

# Arguments:
- `Euv`: the UV energy scale of the spectral density. parameter for DLR grids.
- `rtol`: tolerance absolute error. parameter for DLR grids.
- `param`: parameters of ElectronGas.

# Return
- ``F(\\omega_n, k)=0`` as a `GreenFunc.MeshArray`
- the dynamical part of `R` in the imaginary-time space as a `GreenFunc.MeshArray`, which is set to be zero.
- the instant part of `R` as a `GreenFunc.MeshArray`, ``R_{\\mathrm{ins}}=1``.
"""
# function initFR(Euv, rtol, param)
#     @unpack β = param
#     Ω_c = freq_sep

#     wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
#     R_freq = GreenFunc.MeshArray(wn_mesh; dtype=Float64)
#     F_freq = similar(R_freq)
#     for ind in eachindex(R_freq)
#         ni = ind[1]
#         ωn = wn_mesh[ni]
#         R_freq[ni] = 0.0 #(1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2))
#     end
#     R_imt = real(R_freq |> to_dlr |> to_imtime)
#     #R_ins = GreenFunc.MeshArray([1]; dtype=Float64, data=zeros(1))
#     R_ins = GreenFunc.MeshArray([1]; dtype=Float64, data=ones(1))
#     return F_freq, R_imt, R_ins
# end


"""
    function calcF!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)

Calculation of the Bethe-Slapter amplitude `F` from the product of single-particle Green's function `G2` 
and the dynamical and instant parts of `R`, `R_ins`. Compute in frequency space to avoid \\tau integration.
If F only depends on frequency,
```math
    F(\\omega_n) = G^{(2)}(\\omega_n) [R(\\omega_n)+R_{\\mathrm{ins}}]
```
Otherwise,
```math
    F(\\omega_n, k) = G^{(2)}(\\omega_n, k) [R(\\omega_n, k)+R_{\\mathrm{ins}}(k)]
```
"""
# function calcF!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)
#     R_freq = R |> to_dlr |> to_imfreq
#     # algebraic in frequency space
#     if length(R.mesh) == 1
#         for ind in eachindex(F)
#             F[ind] = real(R_freq[ind] + R_ins[1]) * G2[ind]
#         end
#     else
#         for ind in eachindex(F)
#             F[ind] = real(R_freq[ind] + R_ins[1, ind[2]]) * G2[ind]
#         end
#     end
# end
# function calcF_fromfreq!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)
#     R_freq = R
#     # algebraic in frequency space
#     if length(R.mesh) == 1
#         for ind in eachindex(F)
#             F[ind] = real(R_freq[ind] + R_ins[1]) * G2[ind]
#         end
#     else
#         for ind in eachindex(F)
#             F[ind] = real(R_freq[ind] + R_ins[1, ind[2]]) * G2[ind]
#         end
#     end
# end

function calcF_resum!(F::GreenFunc.MeshArray, F_fs::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray, Πs::GreenFunc.MeshArray; ikF)
    R_freq = R |> to_dlr |> to_imfreq
    # algebraic in frequency space
    for ind in eachindex(F)
        F[ind] = real(R_freq[ind] + R_ins[1, ind[2]]) * G2[ind]
    end

    for ind in eachindex(F_fs)
        F_fs[ind] = -real(R_freq[ind[1], ikF]
                          +
                          R_ins[1, ikF]) * Πs[ind]
    end
    # for ind in 1:length(F_fs.mesh[1])
    #     F_fs[ind, 1] = -real(R_freq[ind, ikF]
    #                          +
    #                          R_ins[1, ikF]) * Πs[ind, 1]
    # end
end

"""
    function calcR!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
        source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite})
 
Calculate ``kR(\\tau, k)`` in three dimensions by given `pF(\\tau, p)` and `kernel`. 
Compute in imaginary time space to aviod frequency convolution.
```math
    kR(\\tau, k) = k\\eta(\\tau, k) - \\int \\frac{dp}{4\\pi^2} H(\\tau,k,p) pF(\\tau,p),
```
where ``H(\\tau,k,p)\\equiv kp W(\\tau,k,p)`` is the helper function of interaction (see [Legendre Decomposition of Interaction](https://numericaleft.github.io/ElectronGas.jl/dev/manual/legendreinteraction/)) 
and the kernel argument. The dynamical `source` ``k\\eta(\\tau, k)`` will be added if it is given as `GreenFunc.MeshArray`.
"""
# function calcR!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
#     source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite})
#     kgrid = F.mesh[2]
#     F_dlr = F |> to_dlr
#     # switch to τ space
#     F_imt = real(F_dlr |> to_imtime)
#     F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)

#     for ind in eachindex(R)
#         # for each τ, k, integrate over q
#         τi, ki = ind[1], ind[2]
#         # interpolate F to q grid of given k
#         Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
#         integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq
#         R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
#         if τi == 1
#             # same for instant part
#             Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
#             integrand = view(kernel_ins, ki, 1:qgrids[ki].size) .* Fq
#             R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
#         end
#     end
#     !(source isa Nothing) && (R += source)
# end

function calcR_resum!(F::GreenFunc.MeshArray, F_fs::GreenFunc.MeshArray,
    R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
    source::Union{Nothing,GreenFunc.MeshArray},
    kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite}; iqFs)
    kgrid = F.mesh[2]
    F_dlr = F |> to_dlr
    # switch to τ space
    F_imt = real(F_dlr |> to_imtime)
    F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)

    Ffs_dlr = F_fs |> to_dlr
    # switch to τ space
    Ffs_imt = real(Ffs_dlr |> to_imtime)
    Ffs_ins = real(dlr_to_imtime(Ffs_dlr, [F_fs.mesh[1].representation.β,])) * (-1)

    for ind in eachindex(R)
        # for each τ, k, integrate over q
        τi, ki = ind[1], ind[2]
        # interpolate F to q grid of given k
        Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
        integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq
        R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        R[τi, ki] += kernel[ki, iqFs[ki], τi] * Ffs_imt[τi, 1] ./ (-4 * π * π)
        if τi == 1
            # same for instant part
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
            integrand = view(kernel_ins, ki, 1:qgrids[ki].size) .* Fq
            R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
            R_ins[1, ki] += kernel_ins[ki, iqFs[ki]] * Ffs_ins[1, 1] ./ (-4 * π * π)
        end
    end
    !(source isa Nothing) && (R += source)
end


"""
    function calcR!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
        source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins)
 
Calculate ``R(\\tau)``  by given `F(\\tau)` and `kernel`. 
Compute in imaginary time space to aviod frequency convolution.
```math
    R(\\tau) = \\eta(\\tau) -  W(\\tau) F(\\tau),
```
where ``W(\\tau)`` is the kernel.
The dynamical `source` ``\\eta(\\tau)`` will be added if it is given as `GreenFunc.MeshArray`.
"""
# function calcR!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
#     source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins)
#     F_dlr = F |> to_dlr
#     # switch to τ space
#     F_imt = real(F_dlr |> to_imtime)
#     F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)

#     R[:] = -F_imt[:] .* kernel[:]
#     #print("$(F_ins[:])\n")
#     R_ins[:] = -F_ins[:] .* kernel_ins[:]
#     #print("$(R_ins[:])\n")
#     !(source isa Nothing) && (R += source)
# end


"""
    function calcR_2d!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,  
        source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite})

Calculate ``kR(\\tau, k)`` in two dimensions by given `pF(p)` and `kernel`
Compute in imaginary time space to aviod frequency convolution.
```math
    kR(\\tau, k) = k\\eta(\\tau, k) - \\int \\frac{pdp}{4\\pi^2} W(\\tau,k,p) pF(\\tau,p),
```
where ``W`` is the interaction and the kernel argument.
The dynamical `source` ``k\\eta(\\tau, k)`` will be added if it is given as `GreenFunc.MeshArray`.
"""
# function calcR_2d!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
#     source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite})
#     # similar to 3d
#     kgrid = F.mesh[2]
#     F_dlr = F |> to_dlr
#     F_imt = real(F_dlr |> to_imtime)
#     F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)
#     for ind in eachindex(R)
#         τi, ki = ind[1], ind[2]
#         k = R.mesh[2][ki]
#         Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
#         integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq .* k
#         R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
#         if τi == 1
#             Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
#             integrand = view(kernel_ins, ki, 1:qgrids[ki].size) .* Fq .* k
#             R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
#         end
#     end
#     !(source isa Nothing) && (R += source)
# end

# function Sigma_γ(n, param)
#     @unpack β, γ, g = param
#     ωn = (2n + 1) * π / β
#     K = (abs(g) * β / 2 / π)^γ
#     A = 0
#     for i in 1:abs(n)
#         A += 2.0 / i^γ
#     end
#     Σ = abs(ωn) #+ π/β*K*A*sign(2n+1)
#     return Σ
# end


"""
    function G02wrapped(Euv, rtol, param)

Returns the inverse self-energy in gamma model.
```math
   \\Sigma^{-1}(\\omega_n) = 1/|\\omega_n+\\pi T K A sgn(2n+1)|
```
where g is the interaction strength, and T is the temperature. Here ``K=\\left( \\frac{g}{2\\pi T} \\right)^{γ}``, ``A=\\sum_{n=1}^{\\infty} \\frac{2}{n^{\\gamma}}``.
"""
# function G2γwrapped(Euv, rtol, param)
#     @unpack β, Ωcut = param
#     wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
#     print("max frequency $(wn_mesh[end]*β)\n")
#     green = GreenFunc.MeshArray(wn_mesh; dtype=Float64)
#     for ind in eachindex(green)
#         ni = ind[1]
#         n = wn_mesh.grid[ni]
#         wn = wn_mesh[ni]
#         # if abs(wn_mesh[ni])>1.0
#         #     green[ind] = 0.0
#         # else
#         #     green[ind] = 1.0/abs(Sigma_γ(n,param))
#         # end
#         green[ind] = 1.0 / Sigma_γ(n, param) / (exp((wn - Ωcut) / 0.0001) + 1) # use fermi distribution to enforce a smooth UV cutoff
#     end
#     return green
# end

"""
    function G02wrapped(Euv, rtol, sgrid, param)

Returns the product of two bare single-particle Green's function.
```math
   G_0^{(2)}(\\omega_n, k) = 1/(\\omega_n^2+\\omega^2)
```
where ``\\omega= k^2/(2m)-\\mu``.
"""
# function G02wrapped(Euv, rtol, sgrid, param)
#     @unpack me, β, μ = param

#     wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
#     green = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
#     for ind in eachindex(green)
#         ni, ki = ind[1], ind[2]
#         ωn = wn_mesh[ni]
#         ω = sgrid.grid[ki]^2 / 2 / me - μ
#         green[ind] = 1 / (ωn^2 + ω^2)
#     end
#     return green
# end

function Πs0wrapped(Euv, rtol, param)
    @unpack me, β, μ, kF, EF = param

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    green = GreenFunc.MeshArray(wn_mesh, [1,]; dtype=Float64)
    for ind in eachindex(green)
        ni, ki = ind[1], ind[2]
        ωn = wn_mesh[ni]
        ω_c = 0.03EF
        ω1 = π * param.T
        # green[ind] = π * me / (kF) / abs(ωn) / (1 + (abs(ωn) / (ω_c))^2) * (1 + (ω1 / ω_c)^2)
        green[ind] = 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    end
    return green
end

"""
    function G2wrapped(Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, param)

Returns the product of two single-particle Green's function from given dynamical and instant parts of self-energy (`Σ` and `Σ_ins`).
```math
   G_0^{(2)}(\\omega_n, k) = 1/\\left([\\omega_n - \\mathrm{Im} \\Sigma(\\omega_n, k)]^2+
   [\\omega+ \\mathrm{Re} \\Sigma(\\omega_n, k) - \\Sigma_{\\mathrm{shift}}]^2\\right)
```
where ``\\omega= k^2/(2m)-\\mu`` and ``\\Sigma_{\\mathrm{shift}}=\\mathrm{Re} \\Sigma(\\omega_0, k_F)``.
"""
# function G2wrapped(Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, param)
#     @unpack me, kF, β, μ = param

#     Σ_freq = Σ |> to_dlr |> to_imfreq
#     green = similar(Σ_freq)

#     w0i_label = locate(Σ_freq.mesh[1], 0)
#     kf_label = locate(Σ_freq.mesh[2], kF)
#     Σ_shift = real(Σ_freq[w0i_label, kf_label] + Σ_ins[1, kf_label])

#     for ind in eachindex(green)
#         ni, ki = ind[1], ind[2]
#         ωn = green.mesh[1][ni]
#         ω = green.mesh[2][ki]^2 / 2 / me - μ
#         ΣR, ΣI = real(Σ_freq[ind] + Σ_ins[1, ki] - Σ_shift), imag(Σ_freq[ind])
#         green[ind] = 1 / ((ωn - ΣI)^2 + (ω + ΣR)^2)
#     end

#     return green
# end

"""
    function BSeq_solver(param, G2::GreenFunc.MeshArray, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite},
        Euv; Ntherm=120, rtol=1e-10, α=0.7, source::Union{Nothing,GreenFunc.MeshArray}=nothing,
        source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1], G2.mesh[2]; dtype=Float64, data=ones(1, G2.mesh[2].size))
    )

Bethe-Slapter equation solver by self-consistent iterations.
```math
    R(\\omega_n, k) = \\eta(\\omega_n, k) - \\frac{1}{\\beta} \\sum_m \\frac{d^dp}{(2\\pi)^d} \\Gamma(\\omega_n,k;\\omega_m,p) G^{(2)}(\\omega_m,p)R(\\omega_m,p) ,
```
where ``\\eta(\\omega_n, k)`` is the sourced term, ``\\Gamma(k,\\omega_n;p,\\omega_m)`` is the particle-particle four-point vertex 
with zero incoming momentum and frequency, and ``G^{(2)}(p,\\omega_m)`` is the product of two single-particle Green's function.

# Arguments
- `param`: parameters of ElectronGas.
- `G2`: product of two single-particle Green's function (::GreenFunc.MeshArray).
- `kernel`: dynamical kernel of the Legendre decomposed effective interaction.
- `kernel_ins`: instant part of the the Legendre decomposed effective interaction.
- `qgrids`: momentum grid of kernel (::Vector{CompositeGrid.Composite}).
- `Euv`: the UV energy scale of the spectral density. 
- `Ntherm`: thermalization step. By defalut, `Ntherm=120`.
- `rtol`: tolerance absolute error. By defalut, `rtol=1e-10`.
- `α`: mixing parameter in the self-consistent iteration. By default, `α=0.7`.
- `source`: dynamical part of sourced term in imaginary-time space. By default, `source=nothing`.
- `source_ins`: instant part of sourced term. By default, `source_ins=1`.

# Return 
- Inverse of low-energy linear response `1/R₀` (``R_0=R(\\omega_0, k_F)``)
- Bethe-Slapter amplitude `F` in imaginary-frequency space
```math
    F(\\omega_m, p) = G^{(2)}(\\omega_m,p)R(\\omega_m,p)
```
- dynamical part of `R` in imaginary-time space
- instant part of `R` in imaginary-time space
"""
# function BSeq_solver(param, G2::GreenFunc.MeshArray, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite},
#     Euv; Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8, source::Union{Nothing,GreenFunc.MeshArray}=nothing,
#     source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1], G2.mesh[2]; dtype=Float64, data=ones(1, G2.mesh[2].size)),
#     verbose=false, Ncheck=5, Nmax=10000)

#     if verbose
#         println("atol=$atol,rtol=$rtol")
#     end

#     @unpack dim, kF = param
#     kgrid = G2.mesh[2]
#     kF_label = locate(kgrid, kF)
#     if !(source isa Nothing)
#         @assert source.mesh[1] isa MeshGrids.ImTime "ImTime is expect for the dim = 1 source."
#         for ni in eachindex(F_freq.mesh[1])
#             source[ni, :] .*= kgrid.grid
#         end
#     end
#     source_ins[1, :] .*= kgrid.grid
#     source0 = source_ins[1, kF_label]

#     # Initalize F and R
#     F_freq, R_imt, R_ins = initFR(Euv, G2.mesh[1].representation.rtol, kgrid, param)

#     R0, R0_sum = 1.0, 0.0
#     dR0, dR0_sum = zeros(Float64, (kgrid.size)), zeros(Float64, (kgrid.size))
#     R_sum = zeros(Float64, (R_imt.mesh[1].representation.size, kgrid.size))

#     lamu, lamu0 = 0.0, 1.0
#     n = 0
#     # self-consistent iteration with mixing α
#     # Note!: all quantites about R, F, and source in the `while` loop are multiplied by momentum k.
#     while (true)
#         n = n + 1
#         # switch between imtime and imfreq to avoid convolution
#         # dlr Fourier transform is much faster than brutal force convolution

#         # calculation from imtime R to imfreq F
#         calcF!(F_freq, R_imt, R_ins, G2)

#         # calculation from imfreq F to imtime R
#         if dim == 3
#             calcR!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
#         elseif dim == 2
#             calcR_2d!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
#         else
#             error("Not implemented for $dim dimensions.")
#         end

#         R_kF = real(dlr_to_imfreq(to_dlr(R_imt), [0])[1, kF_label] + R_ins[1, kF_label])
#         # split the low-energy part R0=R(ω₀,kF) and the remaining instant part dR0 for iterative calcualtions 
#         R0_sum = source0 + R_kF + R0_sum * α
#         dR0_sum = view(R_ins, 1, :) + view(source_ins, 1, :) .- (source0 + R_kF) + dR0_sum .* α
#         R0 = R0_sum * (1 - α)
#         dR0 = dR0_sum .* (1 - α)
#         R_ins[1, :] = dR0 .+ R0

#         # iterative calculation of the dynamical part 
#         R_sum = view(R_imt, :, :) + R_sum .* α
#         R_imt[:, :] = R_sum .* (1 - α)
#         @debug "R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))"
#         # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

#         # record lamu=1/R0 if iterative step n > Ntherm
#         if n > Ntherm && (n % Ncheck == 1)
#             lamu = -1 / (1 + R_kF / kF)
#             if lamu > 0
#                 # this condition does not necessarily mean something wrong
#                 # it only indicates lamu is not converge to correct sign within Ntherm steps
#                 # normally α>0.8 guarantees convergence, then it means Ntherm is too small
#                 @warn ("α = $α or Ntherm=$Ntherm is too small!")
#             end
#             # err = abs(lamu - lamu0)
#             # Exit the loop if the iteration converges
#             isapprox(lamu, lamu0, rtol=rtol, atol=atol) && break
#             n > Nmax && break

#             lamu0 = lamu
#             if verbose
#                 println("lamu=$lamu")
#             end
#         end
#     end
#     println("α = $α, iteration step: $n")
#     lamu = -kF / R0   # calculate 1/R0
#     # calculate real physical quantites F and R
#     for ni in eachindex(F_freq.mesh[1])
#         F_freq[ni, :] ./= kgrid.grid
#         R_imt[ni, :] ./= kgrid.grid
#     end
#     R_ins[1, :] ./= kgrid.grid

#     return lamu, F_freq, R_imt, R_ins
# end

function BSeq_solver_resum(param, G2::GreenFunc.MeshArray, Πs::GreenFunc.MeshArray,
    kernel, kernel_ins,
    qgrids::Vector{CompositeGrid.Composite},
    Euv; Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8,
    source::Union{Nothing,GreenFunc.MeshArray}=nothing,
    source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1], G2.mesh[2];
        dtype=Float64, data=ones(1, G2.mesh[2].size)),
    verbose=false, Ncheck=5, Nmax=10000)

    if verbose
        println("atol=$atol,rtol=$rtol")
    end

    @unpack dim, kF = param
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)
    ikF = kF_label
    iqFs = [locate(qgrids[ki], kF) for ki in 1:kgrid.size]
    println("kF=$kF")
    # println("ikF=$ikF")
    # println("iqFs=$iqFs")

    if !(source isa Nothing)
        @assert source.mesh[1] isa MeshGrids.ImTime "ImTime is expect for the dim = 1 source."
        for ni in eachindex(F_freq.mesh[1])
            source[ni, :] .*= kgrid.grid
        end
    end
    source_ins[1, :] .*= kgrid.grid
    source0 = source_ins[1, kF_label]

    # Initalize F and R
    F_freq, F_fs, R_imt, R_ins = initFR_resum(Euv, G2.mesh[1].representation.rtol, kgrid, param)

    R0, R0_sum = 1.0, 0.0
    dR0, dR0_sum = zeros(Float64, (kgrid.size)), zeros(Float64, (kgrid.size))
    R_sum = zeros(Float64, (R_imt.mesh[1].representation.size, kgrid.size))

    lamu, lamu0 = 0.0, 1.0
    n = 0
    # self-consistent iteration with mixing α
    # Note!: all quantites about R, F, and source in the `while` loop are multiplied by momentum k.
    while (true)
        n = n + 1
        # switch between imtime and imfreq to avoid convolution
        # dlr Fourier transform is much faster than brutal force convolution

        # calculation from imtime R to imfreq F
        calcF_resum!(F_freq, F_fs, R_imt, R_ins, G2, Πs; ikF=ikF)

        # calculation from imfreq F to imtime R
        if dim == 3
            calcR_resum!(F_freq, F_fs, R_imt, R_ins, source, kernel, kernel_ins, qgrids; iqFs=iqFs)
            # elseif dim == 2
            #     calcR_2d!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
        else
            error("Not implemented for $dim dimensions.")
        end

        R_kF = real(dlr_to_imfreq(to_dlr(R_imt), [0])[1, kF_label] + R_ins[1, kF_label])
        # split the low-energy part R0=R(ω₀,kF) and the remaining instant part dR0 for iterative calcualtions 
        R0_sum = source0 + R_kF + R0_sum * α
        dR0_sum = view(R_ins, 1, :) + view(source_ins, 1, :) .- (source0 + R_kF) + dR0_sum .* α
        R0 = R0_sum * (1 - α)
        dR0 = dR0_sum .* (1 - α)
        R_ins[1, :] = dR0 .+ R0

        # iterative calculation of the dynamical part 
        R_sum = view(R_imt, :, :) + R_sum .* α
        R_imt[:, :] = R_sum .* (1 - α)
        @debug "R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))"
        # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

        # record lamu=1/R0 if iterative step n > Ntherm
        if n > Ntherm && (n % Ncheck == 1)
            lamu = -1 / (1 + R_kF / kF)
            if lamu > 0
                # this condition does not necessarily mean something wrong
                # it only indicates lamu is not converge to correct sign within Ntherm steps
                # normally α>0.8 guarantees convergence, then it means Ntherm is too small
                @warn ("α = $α or Ntherm=$Ntherm is too small!")
            end
            # err = abs(lamu - lamu0)
            # Exit the loop if the iteration converges
            isapprox(lamu, lamu0, rtol=rtol, atol=atol) && break
            n > Nmax && break

            lamu0 = lamu
            if verbose
                println("lamu=$lamu")
            end
        end
    end
    println("α = $α, iteration step: $n")
    lamu = -kF / R0   # calculate 1/R0
    # calculate real physical quantites F and R
    for ni in eachindex(F_freq.mesh[1])
        F_freq[ni, :] ./= kgrid.grid
        R_imt[ni, :] ./= kgrid.grid
    end
    R_ins[1, :] ./= kgrid.grid
    F_fs[:, 1] ./ kF

    println("R(∞)=$(R_ins[ikF])")
    return lamu, F_fs, F_freq, R_imt, R_ins
end

"""
    function BSeq_solver(param, G2::GreenFunc.MeshArray, kernel, kernel_ins,
        Euv; Ntherm=120, rtol=1e-10, α=0.7, source::Union{Nothing,GreenFunc.MeshArray}=nothing,
        source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1]; dtype=Float64, data=ones(1))
    )

Bethe-Slapter equation solver by self-consistent iterations (no momentum dependence).
```math
    R(\\omega_n) = \\eta(\\omega_n) - \\frac{1}{\\beta} \\sum_m \\Gamma(\\omega_n;\\omega_m) G^{(2)}(\\omega_m)R(\\omega_m) ,
```
where ``\\eta(\\omega_n)`` is the sourced term, ``\\Gamma(\\omega_n;\\omega_m)`` is the particle-particle four-point vertex 
with zero incoming momentum and frequency, and ``G^{(2)}(\\omega_m)`` is the product of two single-particle Green's function. In γ-model, ``G^{(2)}(\\omega_m)=1/\\Sigma(\\omega_m)``.

# Arguments
- `param`: parameters of ElectronGas.
- `G2`: product of two single-particle Green's function (::GreenFunc.MeshArray).
- `kernel`: dynamical kernel of the Legendre decomposed effective interaction.
- `kernel_ins`: instant part of the the Legendre decomposed effective interaction.
- `Euv`: the UV energy scale of the spectral density. 
- `Ntherm`: thermalization step. By defalut, `Ntherm=120`.
- `rtol`: tolerance absolute error. By defalut, `rtol=1e-10`.
- `α`: mixing parameter in the self-consistent iteration. By default, `α=0.7`.
- `source`: dynamical part of sourced term in imaginary-time space. By default, `source=nothing`.
- `source_ins`: instant part of sourced term. By default, `source_ins=1`.

# Return 
- Inverse of low-energy linear response `1/R₀` (``R_0=R(\\omega_0)``)
- Bethe-Slapter amplitude `F` in imaginary-frequency space
```math
    F(\\omega_m) = G^{(2)}(\\omega_m)R(\\omega_m)
```
- dynamical part of `R` in imaginary-time space
- instant part of `R` in imaginary-time space
"""
# function BSeq_solver(param, G2::GreenFunc.MeshArray, kernel, kernel_ins,
#     Euv; Ntherm=30, rtol=1e-10, atol=1e-10, α=0.8, source::Union{Nothing,GreenFunc.MeshArray}=nothing,
#     source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1]; dtype=Float64, data=ones(1)),
#     verbose=false, Ncheck=5, Nmax=10000,
#     delta_correction=G2 * 0.0) # delta_correction corresponds to additional term in s matrix
#     @unpack β = param
#     if verbose
#         println("atol=$atol,rtol=$rtol")
#     end

#     if !(source isa Nothing)
#         @assert source.mesh[1] isa MeshGrids.ImTime "ImTime is expect for the dim = 1 source."
#     end
#     source0 = source_ins[1]

#     # Initalize F and R
#     # print("G2 $(G2.data)\n $(G2.mesh[1].grid)\n")
#     G_compare = G2 |> to_dlr |> to_imtime
#     # print("G_compare $(G_compare.data)\n")
#     F_freq, R_imt, R_ins = initFR(Euv, G2.mesh[1].representation.rtol, param)
#     R_freq = R_imt |> to_dlr |> to_imfreq

#     #print("$(typeof(R_imt))")
#     R0, R0_sum = 1.0, 0.0
#     dR0, dR0_sum = 0.0, 0.0
#     R_sum = zeros(Float64, (R_freq.mesh[1].representation.size))
#     lamu, lamu0 = 0.0, 1.0
#     n = 0
#     # self-consistent iteration with mixing α
#     # Note!: all quantites about R, F, and source in the `while` loop are multiplied by momentum k.
#     while (true)
#         n = n + 1
#         # switch between imtime and imfreq to avoid convolution
#         # dlr Fourier transform is much faster than brutal force convolution

#         # calculation from imtime R to imfreq F
#         calcF!(F_freq, R_freq, R_ins, G2)

#         # calculation from imfreq F to imtime R
#         calcR!(F_freq, R_imt, R_ins, source, kernel, kernel_ins)
#         R_freq = (R_imt |> to_dlr |> to_imfreq) + (R_freq) .* delta_correction + R_ins[1] * delta_correction
#         # if n == 1
#         #     print("R_ins $(R_ins.data) \n ")
#         # end
#         R_ω0 = real(R_freq[1]) + R_ins[1]
#         # split the low-energy part R0=R(ω₀,kF) and the remaining instant part dR0 for iterative calcualtions 
#         R0_sum = (source0 + R_ω0) + R0_sum * α
#         dR0_sum = view(R_ins, 1) + view(source_ins, 1) .- (source0 + R_ω0) + dR0_sum .* α
#         R0 = R0_sum * (1 - α)
#         dR0 = dR0_sum .* (1 - α)
#         R_ins[1] = R0 + dR0

#         # iterative calculation of the dynamical part 
#         R_sum = view(R_freq, :) + R_sum .* α
#         R_freq[:] = R_sum .* (1 - α)
#         @debug "R(ω0) = $R_ω0, 1/R0 = $(-1/R0)  ($(-1/(1 + R_ω0)))"
#         # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

#         # record lamu=1/R0 if iterative step n > Ntherm
#         if n > Ntherm && (n % Ncheck == 1)
#             lamu = -1 / (1 + R_ω0)
#             if lamu > 0
#                 # this condition does not necessarily mean something wrong
#                 # it only indicates lamu is not converge to correct sign within Ntherm steps
#                 # normally α>0.8 guarantees convergence, then it means Ntherm is too small
#                 @warn ("α = $α or Ntherm=$Ntherm is too small!")
#             end
#             # err = abs(lamu - lamu0)
#             # Exit the loop if the iteration converges
#             isapprox(lamu, lamu0, rtol=rtol, atol=atol) && break
#             n > Nmax && println("[WARNING] reach Nmax=$Nmax") && break

#             lamu0 = lamu
#             if verbose
#                 println("lamu=$lamu")
#             end
#         end
#     end
#     println("α = $α, iteration step: $n")
#     lamu = 1.0 / R0   # calculate 1/R0
#     # calculate real physical quantites F and R

#     return lamu, F_freq, R_imt, R_ins
# end




"""
    function linearResponse(param, channel::Int; Euv=100 * param.EF, rtol=1e-10,
        maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8, sigmatype=:none, int_type=:rpa, α=0.7)

Implmentation of Cooper-pair linear response approach.

# Arguments:
- `param`: parameters of ElectronGas.
- `channel::Int`: orbital angular channel (0: s-wave, 1: p-wave, ...)
- `Euv`: the UV energy scale of the spectral density. By default, `Euv=100EF`.
- `rtol`: tolerance absolute error. By defalut, `rtol=1e-10`.
- `maxK`: maximum momentum of kgrid and qgrids. By default, `maxK=10kF`.
- `minK`: minimum interval of panel kgrid and qgrids. By default, `minK=1e-7kF`.
- `Nk`: number of grid points of panel kgrid and qgrids. By defalut, `Nk=8`.
- `order`: number of grid points of subgrid of kgrid and qgrids. By defalut, `order=8`.
- `sigmatype`: type of fermionic self-energy. (no self-energy :none, G0W0 approximation :g0w0)
- `int_type`: type of effective interaction. By default, `int_type=:rpa`.
- `α`: mixing parameter in the self-consistent iteration. By default, `α=0.7`.

# Return:
- Inverse of low-energy linear response `1/R₀` (``R_0=R(\\omega_0, k_F)``)
- Linear response `R(ωₙ, k)`, which is calculated by the Bethe-Slapter-type equation
```math
    R(\\omega_n, k) = 1 - \\frac{1}{\\beta} \\sum_m \\int \\frac{d^dp}{(2\\pi)^d} \\Gamma(\\omega_n,k;\\omega_m,p) G^{(2)}(\\omega_m,p)R(\\omega_m,p)
```
where ``1`` is the default sourced term, ``\\Gamma(k,\\omega_n;p,\\omega_m)`` is the particle-particle four-point vertex 
with zero incoming momentum and frequency, and ``G^{(2)}(p,\\omega_m)`` is the product of two single-particle Green's function.
"""
function linearResponse(param, channel::Int;
    Euv=100 * param.EF, rtol=1e-10, atol=1e-10,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8,
    Vph::Union{Function,Nothing}=nothing, sigmatype=:none, int_type=:rpa,
    α=0.8, verbose=false, Ntherm=30, Nmax=10000,
    resum=false,
    issave=false, uid=1, dir="./", kwargs...)
    println("resum=$resum")
    @unpack dim, rs, β, kF = param
    if verbose
        println("atol=$atol,rtol=$rtol")
    end
    # prepare Legendre decomposed effective interaction
    if dim == 3
        if channel == 0
            # nmax = G2.mesh[1].grid[end]
            # delta_correction = Δω_correction(G2, wsph, a2f_iso, nmax)
            @time W = LegendreInteraction.DCKernel0(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        else
            @time W = LegendreInteraction.DCKernel_old(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph, kwargs...)
        end
    elseif dim == 2
        @time W = LegendreInteraction.DCKernel_2d(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
            minK=minK, order=order, int_type=int_type, channel=channel, Vph=Vph)
    else
        error("No support for $dim dimension!")
    end

    fdlr = Lehmann.DLRGrid(Euv, β, rtol, true, :pha)
    bdlr = W.dlrGrid
    kgrid = W.kgrid
    qgrids = W.qgrids

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_ins = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3))

    kF_label = locate(kgrid, kF)
    qF_label = locate(qgrids[kF_label], kF)
    println("static kernel at (kF, kF):$(kernel_ins[kF_label, qF_label])")

    # prepare G2 as sigmatype
    if sigmatype == :none
        G2 = BSeq.G02wrapped(Euv, rtol, kgrid, param)
    elseif sigmatype == :g0w0
        @unpack me, β, μ = param
        wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
        Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        # self energy should be converted to proper frequency grid
        Σ_dlr = Σ |> to_dlr
        Σ_wn = dlr_to_imfreq(Σ_dlr, wn_mesh)
        G2 = BSeq.G2wrapped(Σ_wn, Σ_ins, param)
    end

    # calculate F, R by Bethe-Slapter iteration.
    if resum
        Πs = Πs0wrapped(Euv, rtol, param)
        lamu, F_fs, F_freq, R_imt, R_ins = BSeq_solver_resum(param, G2, Πs, kernel, kernel_ins, qgrids, Euv;
            rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
    else
        lamu, F_freq, R_imt, R_ins = BSeq.BSeq_solver(param, G2, kernel, kernel_ins, qgrids, Euv;
            rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
    end
    println("1/R₀ = $lamu")

    R_freq = R_imt |> to_dlr |> to_imfreq
    if issave
        fname = "PCFdata_$(uid).jld2"
        jldopen(dir * fname, "w") do file
            file["param"] = param
            file["lamu"] = lamu
            file["F_freq"] = F_freq
            file["R_ins"] = R_ins
            file["R_freq"] = R_freq
        end
    end
    # println(view(R_freq, :, kF_label))
    return lamu, R_freq, F_freq, R_ins
end


"""
    function linearResponse(param; Euv=100 * param.EF, rtol=1e-10,
        maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8, α=0.8)

Implmentation of γ-model linear response approach.

# Arguments:
- `param`: parameters of ElectronGas.
- `Euv`: the UV energy scale of the spectral density. By default, `Euv=100EF`.
- `rtol`: tolerance absolute error. By defalut, `rtol=1e-10`.
- `maxK`: maximum momentum of kgrid and qgrids. By default, `maxK=10kF`.
- `minK`: minimum interval of panel kgrid and qgrids. By default, `minK=1e-7kF`.
- `Nk`: number of grid points of panel kgrid and qgrids. By defalut, `Nk=8`.
- `order`: number of grid points of subgrid of kgrid and qgrids. By defalut, `order=8`.
- `α`: mixing parameter in the self-consistent iteration. By default, `α=0.7`.

# Return:
- Inverse of low-energy linear response `1/R₀` (``R_0=R(\\omega_0)``)
- Linear response `R(ωₙ)`, which is calculated by the Bethe-Slapter-type equation
```math
    R(\\omega_n) = 1 - \\frac{1}{\\beta} \\sum_m\\Gamma(\\omega_n;\\omega_m) G^{(2)}(\\omega_m)R(\\omega_m)
```
where ``1`` is the default sourced term, ``\\Gamma(\\omega_n;\\omega_m)`` is the particle-particle four-point vertex 
with zero incoming momentum and frequency, and ``G^{(2)}(\\omega_m)`` is the product of two single-particle Green's function.
"""
# function linearResponse(param; Euv=100 * param.EF, rtol=1e-10, atol=1e-10, α=0.7, verbose=false, Ntherm=30, Nmax=10000)
#     @unpack g, β, γ = param
#     if verbose
#         println("atol=$atol,rtol=$rtol")
#     end
#     # prepare Legendre decomposed effective interaction
#     @time kernel_freq, kernel_ins = Interaction.gamma_wrapped(Euv, rtol, param)
#     G2 = G2γwrapped(Euv, rtol, param)
#     fdlr = Lehmann.DLRGrid(Euv, β, rtol, FERMION, :pha)
#     # bdlr = kernel_freq.mesh[1].dlrGrid
#     @assert G2.mesh[1].grid == fdlr.n
#     # prepare kernel, interpolate into τ-space with fdlr.τ
#     # kernel =  kernel_freq |> to_dlr |>to_imtime
#     kernel = dlr_to_imtime(to_dlr(kernel_freq), fdlr.τ)

#     # calculate F, R by Bethe-Slapter iteration.
#     lamu, F_freq, R_imt, R_ins = BSeq_solver(param, G2, kernel, kernel_ins, Euv;
#         rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, Nmax=Nmax)
#     println("1/R₀ = $lamu")

#     R_freq = R_imt |> to_dlr |> to_imfreq
#     # println(view(R_freq, :, kF_label))
#     return lamu, R_freq, F_freq
# end

# include("gap_function.jl")

end

