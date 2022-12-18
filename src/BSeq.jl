"""
Bethe-Slapter-type equation solver and the application to Cooper-pair linear response approach. 
"""
module BSeq

using ..Parameter, ..Convention, ..LegendreInteraction
using ..Parameters, ..GreenFunc, ..Lehmann, ..CompositeGrids

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
function initFR(Euv, rtol, sgrid, param)
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
    return F_freq, R_imt, R_ins
end

"""
    function calcF!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)

Calculation of the Bethe-Slapter amplitude `F` from the product of single-particle Green's function `G2` 
and the dynamical and instant parts of `R`, `R_ins`. Compute in frequency space to avoid \\tau integration.
```math
    F(\\omega_n, k) = G^{(2)}(\\omega_n, k) [R(\\omega_n, k)+R_{\\mathrm{ins}}(k)]
```
"""
function calcF!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)
    R_freq = R |> to_dlr |> to_imfreq
    # algebraic in frequency space
    for ind in eachindex(F)
        F[ind] = real(R_freq[ind] + R_ins[1, ind[2]]) * G2[ind]
    end
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
function calcR!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
    source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite})
    kgrid = F.mesh[2]
    F_dlr = F |> to_dlr
    # switch to τ space
    F_imt = real(F_dlr |> to_imtime)
    F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)

    for ind in eachindex(R)
        # for each τ, k, integrate over q
        τi, ki = ind[1], ind[2]
        # interpolate F to q grid of given k
        Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
        integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq
        R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        if τi == 1
            # same for instant part
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
            integrand = view(kernel_ins, ki, 1:qgrids[ki].size) .* Fq
            R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        end
    end
    !(source isa Nothing) && (R += source)
end

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
function calcR_2d!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray,
    source::Union{Nothing,GreenFunc.MeshArray}, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite})
    # similar to 3d
    kgrid = F.mesh[2]
    F_dlr = F |> to_dlr
    F_imt = real(F_dlr |> to_imtime)
    F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)
    for ind in eachindex(R)
        τi, ki = ind[1], ind[2]
        k = R.mesh[2][ki]
        Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
        integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq .* k
        R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        if τi == 1
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
            integrand = view(kernel_ins, ki, 1:qgrids[ki].size) .* Fq .* k
            R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        end
    end
    !(source isa Nothing) && (R += source)
end

"""
    function G02wrapped(Euv, rtol, sgrid, param)

Returns the product of two bare single-particle Green's function.
```math
   G_0^{(2)}(\\omega_n, k) = 1/(\\omega_n^2+\\omega^2)
```
where ``\\omega= k^2/(2m)-\\mu``.
"""
function G02wrapped(Euv, rtol, sgrid, param)
    @unpack me, β, μ = param

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    green = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
    for ind in eachindex(green)
        ni, ki = ind[1], ind[2]
        ωn = wn_mesh[ni]
        ω = sgrid.grid[ki]^2 / 2 / me - μ
        green[ind] = 1 / (ωn^2 + ω^2)
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
function G2wrapped(Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, param)
    @unpack me, kF, β, μ = param

    Σ_freq = Σ |> to_dlr |> to_imfreq
    green = similar(Σ_freq)

    w0i_label = locate(Σ_freq.mesh[1], 0)
    kf_label = locate(Σ_freq.mesh[2], kF)
    Σ_shift = real(Σ_freq[w0i_label, kf_label] + Σ_ins[1, kf_label])

    for ind in eachindex(green)
        ni, ki = ind[1], ind[2]
        ωn = green.mesh[1][ni]
        ω = green.mesh[2][ki]^2 / 2 / me - μ
        ΣR, ΣI = real(Σ_freq[ind] + Σ_ins[1, ki] - Σ_shift), imag(Σ_freq[ind])
        green[ind] = 1 / ((ωn - ΣI)^2 + (ω + ΣR)^2)
    end

    return green
end

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
function BSeq_solver(param, G2::GreenFunc.MeshArray, kernel, kernel_ins, qgrids::Vector{CompositeGrid.Composite},
    Euv; Ntherm=120, rtol=1e-10, α=0.7, source::Union{Nothing,GreenFunc.MeshArray}=nothing,
    source_ins::GreenFunc.MeshArray=GreenFunc.MeshArray([1], G2.mesh[2]; dtype=Float64, data=ones(1, G2.mesh[2].size))
)
    @unpack dim, kF = param
    kgrid = G2.mesh[2]
    kF_label = locate(kgrid, kF)
    if !(source isa Nothing)
        @assert source.mesh[1] isa MeshGrids.ImTime "ImTime is expect for the dim = 1 source."
        for ni in eachindex(F_freq.mesh[1])
            source[ni, :] .*= kgrid.grid
        end
    end
    source_ins[1, :] .*= kgrid.grid
    source0 = source_ins[1, kF_label]

    # Initalize F and R
    F_freq, R_imt, R_ins = initFR(Euv, G2.mesh[1].representation.rtol, kgrid, param)

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
        calcF!(F_freq, R_imt, R_ins, G2)

        # calculation from imfreq F to imtime R
        if dim == 3
            calcR!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
        elseif dim == 2
            calcR_2d!(F_freq, R_imt, R_ins, source, kernel, kernel_ins, qgrids)
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
        if n > Ntherm
            lamu = -1 / (1 + R_kF / kF)
            lamu > 0 && error("α = $α is too small!")
            err = abs(lamu - lamu0) / abs(lamu + EPS)
            # Exit the loop if the iteraction converges
            lamu >= lamu0 > -1 && err < rtol && break
            lamu0 <= lamu < -1 && err < rtol && break

            lamu0 = lamu
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

    return lamu, F_freq, R_imt, R_ins
end

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
function linearResponse(param, channel::Int; Euv=100 * param.EF, rtol=1e-10,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8, sigmatype=:none, int_type=:rpa, α=0.7)
    @unpack dim, rs, β, kF = param

    # prepare Legendre decomposed effective interaction
    if dim == 3
        if channel == 0
            @time W = LegendreInteraction.DCKernel0(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel)
        else
            @time W = LegendreInteraction.DCKernel_old(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
                minK=minK, order=order, int_type=int_type, channel=channel)
        end
    elseif dim == 2
        @time W = LegendreInteraction.DCKernel_2d(param; Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK,
            minK=minK, order=order, int_type=int_type, channel=channel)
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
        G2 = G02wrapped(Euv, rtol, kgrid, param)
    elseif sigmatype == :g0w0
        Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        G2 = G2wrapped(Σ, Σ_ins, param)
    end

    # calculate F, R by Bethe-Slapter iteration.
    lamu, F_freq, R_imt, R_ins = BSeq_solver(param, G2, kernel, kernel_ins, qgrids, Euv; rtol=rtol, α=α)
    println("1/R₀ = $lamu")

    R_freq = R_imt |> to_dlr |> to_imfreq
    # println(view(R_freq, :, kF_label))
    return lamu, R_freq
end

end