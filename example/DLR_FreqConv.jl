module FreqConv

using Lehmann, StaticArrays

function fullMat(α, γ, n, β)
    ω = π * (2n + 1) / β

    if α == 0
        return Spectral.kernelΩ(Val(false), Val(:ph), 0, α, β) * Spectral.kernelΩ(Val(true), Val(:pha), n, γ, β) / β
    end
    factor = 4 * α * γ * (-expm1(-α * β)) * (1 + ℯ^(-γ * β))
    f_term = 1.0 / ((α^2 - γ^2 + ω^2)^2 + 4 * γ^2 * ω^2) * (α^2 - γ^2 + ω^2) / (2γ) / (1 + ℯ^(-γ * β)) * (-expm1(-γ * β))
    b_term = 1.0 / ((-α^2 + γ^2 + ω^2)^2 + 4 * α^2 * ω^2) * (-α^2 + γ^2 + ω^2) / (2α) / (-expm1(-α * β)) * (1 + ℯ^(-α * β))
    #    println(factor, "\t", f_term, "\t", b_term)
    return factor * (f_term + b_term)
end

function sumMat(α, γ, n, β, n_c)
    result = 0.0
    for m in -n_c-1:n_c
        result += Spectral.kernelΩ(Val(true), Val(:pha), m, γ, β) * Spectral.kernelΩ(Val(false), Val(:ph), n - m, α, β)
    end

    return result / β
end



struct ConvMat
    bdlr::DLRGrid
    fdlr::DLRGrid

    n_c::Int

    full_mat::Array{Float64,3}
    sep_mat::Array{Float64,3}
    high_mat::Array{Float64,3}

    asw_full::Array{Float64,1}
    asw_low::Array{Float64,1}
    asw_high::Array{Float64,1}

    function ConvMat(bdlr::DLRGrid, fdlr::DLRGrid, n_c::Int)
        sep_mat = zeros(Float64, (fdlr.size, bdlr.size, fdlr.size))
        full_mat = zeros(Float64, (fdlr.size, bdlr.size, fdlr.size))
        β = fdlr.β

        for (ni, n) in enumerate(fdlr.n)
            for (ωi, ω) in enumerate(bdlr.ω)
                for (ξi, ξ) in enumerate(fdlr.ω)
                    # for m in 1:n_c
                    #     sep_mat[ni, ωi, ξi] += Spectral.kernelΩ(Val(true), Val(:pha), m, ξ, β)*Spectral.kernelΩ(Val(false), Val(:ph),n-m, ω, β)
                    # end
                    sep_mat[ni, ωi, ξi] = sumMat(ω, ξ, n, β, n_c)
                    full_mat[ni, ωi, ξi] = fullMat(ω, ξ, n, β)
                end
            end
        end
        high_mat = full_mat .- sep_mat

        asw_full = zeros(Float64, fdlr.size)
        asw_low = zeros(Float64, fdlr.size)
        asw_high = zeros(Float64, fdlr.size)
        for (ξi, ξ) in enumerate(fdlr.ω)
            for m in -n_c-1:n_c
                asw_low[ξi] += Spectral.kernelΩ(Val(true), Val(:pha), m, ξ, β) / β
            end
            asw_full[ξi] = (1 + exp(-ξ * β)) * tanh(ξ * β / 2.0)
        end
        asw_high = asw_full .- asw_low

        return new(bdlr, fdlr, n_c, full_mat, sep_mat, high_mat, asw_full, asw_low, asw_high)
    end

end

function freq_conv(Γ::Vector, F::Vector, CM::ConvMat, type=:full)
    @assert length(Γ) == CM.bdlr.size
    @assert length(F) == CM.fdlr.size

    conv_mat = CM.full_mat
    if type == :low
        conv_mat = CM.sep_mat
    elseif type == :high
        conv_mat -= CM.sep_mat
    end

    Δ = zeros(ComplexF64, CM.fdlr.size)

    γ = matfreq2dlr(CM.bdlr, Γ)
    f = matfreq2dlr(CM.fdlr, F)

    # println(CM.fdlr.n)
    # println(CM.bdlr.ω)
    # println(CM.fdlr.ω)
    for (ni, n) in enumerate(CM.fdlr.n)
        for (ωi, ω) in enumerate(CM.bdlr.ω)
            for (ξi, ξ) in enumerate(CM.fdlr.ω)
                Δ[ni] += γ[ωi] * f[ξi] * conv_mat[ni, ωi, ξi]
            end
        end
    end

    return real(Δ)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    const E_c = π
    const Ω_c = 1.0
    const Ω_c2 = π
    const g = 2.0

    const fEUV = 10
    const bEUV = 10

    const β = 200.0

    function Phonon(ω)
        #return g .* ω .^ 2 ./ (ω .^ 2 .+ Ω_c^2 ) #.- g .* ω .^ 2 ./ (ω .^ 2 .+ Ω_c2^2 )
        return g ./ (ω .^ 2 .+ Ω_c^2)
    end

    function AcorrSepWeight(fdlr, Ω_c)
        β = fdlr.β
        n_c = Base.floor(Int, Ω_c * β / (2π) - 0.5)
        ASW = zeros(Float64, fdlr.size)
        for (ξi, ξ) in enumerate(fdlr.ω)
            for m in 1:n_c
                ASW[ξi] += FreqConv.Spectral.kernelΩ(Val(true), Val(:pha), m, ξ, β)
            end
        end

        return ASW
    end

    fdlr = FreqConv.DLRGrid(fEUV, β, 1e-12, true, :pha)
    bdlr = FreqConv.DLRGrid(bEUV, β, 1e-12, false, :ph)

    println("α:", bdlr.ω)
    println("γ:", fdlr.ω)

    println("ASW:", AcorrSepWeight(fdlr, 100000.0))
    println("ASW_inf:",)
    #    α, ξ = bdlr.ω[2],fdlr.ω[3]
    α, ξ = 0.00073, 0.056
    n = 2
    m = 5
    println("sumMat:", FreqConv.sumMat(α, ξ, n, β, 100), "\t", FreqConv.sumMat(α, ξ, n, β, 1000), "\t", FreqConv.sumMat(α, ξ, n, β, 10000), "\t")
    println("fullMat:", FreqConv.fullMat(α, ξ, n, β))
    println(FreqConv.Spectral.kernelΩ(Val(false), Val(:ph), n - m, α, β) * FreqConv.Spectral.kernelΩ(Val(true), Val(:pha), m, ξ, β))

    #    @assert 1==2 "break"

    n_max = Base.floor(Int, E_c / (2π) * β)
    n_c = Base.floor(Int, Ω_c / (2π) * β)
    println("n_c=$(n_c)")
    N = 33

    cm = FreqConv.ConvMat(bdlr, fdlr, n_c)

    Γ = Phonon(bdlr.ωn)
    ωf = 2π .* (fdlr.n .+ 0.5) ./ β
    F = (1.0 .- 0.1 .* Phonon(ωf)) ./ ωf
    println(Γ)
    println(F)
    γ = FreqConv.matfreq2dlr(bdlr, Γ)
    f = FreqConv.matfreq2dlr(fdlr, F)
    Γt = FreqConv.dlr2tau(bdlr, γ, bdlr.τ)
    Ft = FreqConv.dlr2tau(fdlr, f, fdlr.τ)
    # Γt = FreqConv.MultiPole(:corr, bdlr.τ, β, bEUV)[1]
    # Ft = FreqConv.MultiPole(:acorr, fdlr.τ, β, fEUV)[1]
    # γ = FreqConv.tau2dlr(bdlr, Γt)
    # f = FreqConv.tau2dlr(fdlr, Ft)
    # Γ = FreqConv.dlr2matfreq(bdlr, γ, bdlr.n)
    # F = FreqConv.dlr2matfreq(fdlr, f, fdlr.n)


    Δ = FreqConv.freq_conv(Γ, F, cm, :low)
    println(fdlr.n[1:N])
    println("Δ:", Δ[1:N])

    Δ_comp = zeros(Float64, fdlr.size)
    Γ_mat = zeros(Float64, (fdlr.size, 2 * n_c + 2))
    F_v = zeros(Float64, 2 * n_c + 2)
    # Γ_mat = zeros(Float64, (fdlr.size, n_c+1))
    # F_v = zeros(Float64, n_c+1)

    for (ni, n) in enumerate(fdlr.n)
        Γ_mat[ni, :] = FreqConv.dlr2matfreq(bdlr, γ, [n - m for m in -n_c-1:n_c])
        #Γ_mat[ni, :] = FreqConv.dlr2matfreq(bdlr, γ, [n-m for m in 0:n_c])
    end
    F_v = FreqConv.dlr2matfreq(fdlr, f, [m for m in -n_c-1:n_c])
    #F_v = FreqConv.dlr2matfreq(fdlr, f, [m for m in 0:n_c])

    println(Γ_mat[1, :])
    println(F_v)

    for (ni, n) in enumerate(fdlr.n)
        for m in 1:n_c*2+2
            #for m in 1:n_c*1+1
            # Δ_comp[ni] += test_spec(2π*(n-m)/β) * test_spec(2π*(m+0.5)/β)
            Δ_comp[ni] += Γ_mat[ni, m] * F_v[m] / β
        end
    end

    println("Δ_comp:", Δ_comp[1:N])

    cm = FreqConv.ConvMat(bdlr, fdlr, 1000)

    Δ = FreqConv.freq_conv(Γ, F, cm, :full)
    println(Δ)
    Γft = FreqConv.dlr2tau(bdlr, γ, fdlr.τ)
    Δ2t = Γft .* Ft
    Δ2 = FreqConv.tau2matfreq(fdlr, Δ2t, fdlr.n)
    println(real(Δ2))

    F0 = FreqConv.dlr2tau(fdlr, f, [0.0,])[1]
    F1 = 0.0
    F2 = 0.0
    F3 = 0.0
    for ni in 1:fdlr.size
        global F1 = F1 + f[ni] * cm.asw_low[ni]
        global F2 = F2 + f[ni] * cm.asw_high[ni]
        global F3 = F3 + f[ni] * cm.asw_full[ni]
    end
    println("$(F0), $(F1), $(F2), $(F3)")

end