"""
Solve Bethe-Slapter equation
"""

module BSeq

using ..Parameter, ..Convention, ..LegendreInteraction
using ..Parameters, ..GreenFunc, ..Lehmann, ..CompositeGrids

const freq_sep = 0.01

function initFR(Euv, rtol, sgrid, param)
    @unpack β, me, μ = param
    Ω_c = freq_sep

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    R_freq = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
    F_freq = similar(R_freq)
    for ind in eachindex(R_freq)
        ni, ki = ind[1], ind[2]
        ωn, k = matfreq(wn_mesh)[ni], sgrid[ki]
        e = k^2 / 2 / me - μ
        R_freq[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
    end
    R_imt = real(R_freq |> to_dlr |> to_imtime)
    R_ins = GreenFunc.MeshArray([1], sgrid; dtype=Float64, data=zeros(1, sgrid.size))
    return F_freq, R_imt, R_ins
end

function calcF!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, G2::GreenFunc.MeshArray)
    R_freq = R |> to_dlr |> to_imfreq
    for ind in eachindex(F)
        F[ind] = real(R_freq[ind] + R_ins[1, ind[2]]) * G2[ind]
    end
end

function calcR!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, kernel, kernel_bare, qgrids)
    kgrid = F.mesh[2]
    F_dlr = F |> to_dlr
    F_imt = real(F_dlr |> to_imtime)
    F_ins = real(dlr_to_imtime(F_dlr, [F.mesh[1].representation.β,])) * (-1)
    for ind in eachindex(R)
        τi, ki = ind[1], ind[2]
        Fq = CompositeGrids.Interp.interp1DGrid(view(F_imt, τi, :), kgrid, qgrids[ki].grid)
        integrand = view(kernel, ki, 1:qgrids[ki].size, τi) .* Fq
        R[τi, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        if τi == 1
            Fq = CompositeGrids.Interp.interp1DGrid(view(F_ins, 1, :), kgrid, qgrids[ki].grid)
            integrand = view(kernel_bare, ki, 1:qgrids[ki].size) .* Fq
            R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        end
    end
end

function calcR_2d!(F::GreenFunc.MeshArray, R::GreenFunc.MeshArray, R_ins::GreenFunc.MeshArray, kernel, kernel_bare, qgrids)
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
            integrand = view(kernel_bare, ki, 1:qgrids[ki].size) .* Fq .* k
            R_ins[1, ki] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki]) ./ (-4 * π * π)
        end
    end
end

function G02wrapped(Euv, rtol, sgrid, param)
    @unpack me, β, μ = param

    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    green = GreenFunc.MeshArray(wn_mesh, sgrid; dtype=Float64)
    for ind in eachindex(green)
        ni, ki = ind[1], ind[2]
        ωn = matfreq(wn_mesh)[ni]
        ω = sgrid.grid[ki]^2 / 2 / me - μ
        green[ind] = 1 / (ωn^2 + ω^2)
    end
    return green
end

function G2wrapped(Σ::GreenFunc.MeshArray, Σ_ins::GreenFunc.MeshArray, param)
    @unpack me, kF, β, μ = param

    Σ_freq = Σ |> to_dlr |> to_imfreq
    green = similar(Σ_freq)

    w0i_label = locate(Σ_freq.mesh[1], 0)
    kf_label = locate(Σ_freq.mesh[2], kF)
    Σ_shift = real(Σ_freq[w0i_label, kf_label] + Σ_ins[1, kf_label])

    for ind in eachindex(green)
        ni, ki = ind[1], ind[2]
        ωn = matfreq(green.mesh[1])[ni]
        ω = green.mesh[2][ki]^2 / 2 / me - μ
        ΣR, ΣI = real(Σ_freq[ind] + Σ_ins[1, ki] - Σ_shift), imag(Σ_freq[ind])
        green[ind] = 1 / ((ωn - ΣI)^2 + (ω + ΣR)^2)
    end

    return green
end

function BSeq_solver(param, G2::GreenFunc.MeshArray, kernel, kernel_bare, qgrids, Euv;
    Ntherm=120, rtol=1e-10, α=0.7, source=1, source_kF=1)
    @unpack dim, kF = param
    kgrid = G2.mesh[2]

    F_freq, R_imt, R_ins = initFR(Euv, G2.mesh[1].representation.rtol, kgrid, param)

    kF_label = locate(kgrid, kF)
    R0, R0_sum = 1.0, 0.0
    dR0, dR0_sum = zeros(Float64, (kgrid.size)), zeros(Float64, (kgrid.size))
    R_sum = zeros(Float64, (R_imt.mesh[1].representation.size, kgrid.size))

    lamu, lamu0 = 0.0, 1.0
    n = 0
    while (true)
        n = n + 1

        calcF!(F_freq, R_imt, R_ins, G2)
        if dim == 3
            calcR!(F_freq, R_imt, R_ins, kernel, kernel_bare, qgrids)
        elseif dim == 2
            calcR_2d!(F_freq, R_imt, R_ins, kernel, kernel_bare, qgrids)
        else
            error("Not implemented for $dim dimensions.")
        end

        R_kF = real(dlr_to_imfreq(to_dlr(R_imt), [0])[1, kF_label] + R_ins[1, kF_label])

        R0_sum = kF * source_kF + R_kF + R0_sum * α
        dR0_sum = view(R_ins, 1, :) + kgrid.grid .- (kF * source_kF + R_kF) + dR0_sum .* α
        R0 = R0_sum * (1 - α)
        dR0 = dR0_sum .* (1 - α)
        R_ins[1, :] = dR0 .+ R0

        R_sum = view(R_imt, :, :) + R_sum .* α
        R_imt[:, :] = R_sum .* (1 - α)
        @debug "R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))"
        # println("R(ω0, kF) = $R_kF, 1/R0 = $(-1/R0)  ($(-1/(kF + R_kF)))")

        if n > Ntherm
            lamu = -1 / (1 + R_kF / kF)
            lamu > 0 && error("α = $α is too small!")

            err = abs(lamu - lamu0) / abs(lamu + EPS)
            lamu >= lamu0 > -1 && err < rtol && break
            lamu0 <= lamu < -1 && err < rtol && break
            lamu0 = lamu
        end
    end
    println("α = $α, iteration step: $n")
    lamu = -kF / R0
    return lamu, F_freq, R_imt, R_ins
end

function linearResponse(param, channel::Int; Euv=100 * param.EF, rtol=1e-10,
    maxK=10param.kF, minK=1e-7param.kF, Nk=8, order=8, sigmatype=:none, int_type=:rpa, α=0.7)
    @unpack dim, rs, β, kF = param
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

    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis=3))

    kF_label = locate(kgrid, kF)
    qF_label = locate(qgrids[kF_label], kF)
    println("static kernel at (kF, kF):$(kernel_bare[kF_label, qF_label])")

    if sigmatype == :none
        G2 = G02wrapped(Euv, rtol, kgrid, param)
    elseif dim == 3
        Σ, Σ_ins = SelfEnergy.G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
        G2 = G2wrapped(Σ, Σ_ins, param)
    end

    lamu, F_freq, R_imt, R_ins = BSeq_solver(param, G2, kernel, kernel_bare, qgrids, Euv; rtol=rtol, α=α)
    println("1/R₀ = $lamu")

    R_freq = R_imt |> to_dlr |> to_imfreq
    # println(view(R_freq, :, kF_label))
    return lamu, R_freq
end

end