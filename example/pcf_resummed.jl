using ElectronGas
using JLD2
using ElectronGas.CompositeGrids
using ElectronGas.Parameters
using ElectronGas.GreenFunc

function load_AB(fname)
    f = jldopen(fname, "r")
    param = f["param"]
    A = f["A"]
    B = f["B"]
    return param, A, B
end

function load_B(fname)
    f = jldopen(fname, "r")
    param = f["param"]
    B = f["B"]
    return param, B
end

function load_A(fname)
    f = jldopen(fname, "r")
    param = f["param"]
    A = f["A"]
    return param, A
end

function extend_AB(A, B, param;
    Nk=8, minterval=0.001, order=4)
    # extend A and B to 0
    oldwgrid = A.mesh[1]
    Euv = oldwgrid[end]
    wgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [0.0, Euv], [0.0, 0.1param.EF], Nk, minterval, order)
    newA = GreenFunc.MeshArray(wgrid; dtype=Float64)
    newB = GreenFunc.MeshArray(wgrid, wgrid; dtype=Float64)
    for i in 1:length(wgrid)
        w = wgrid[i]
        if w < oldwgrid[1]
            newA.data[i] = A.data[1]
        else
            newA.data[i] = CompositeGrids.Interp.interp1D(A.data, A.mesh[1], w)
        end
        for j in 1:length(wgrid)
            v = wgrid[j]
            if w < oldwgrid[1] && v < oldwgrid[1]
                newB.data[i, j] = B.data[1, 1]
            elseif w < oldwgrid[1]
                newB.data[i, j] = CompositeGrids.Interp.interp1D(view(B.data, 1, :), B.mesh[2], v)
            elseif v < oldwgrid[1]
                newB.data[i, j] = CompositeGrids.Interp.interp1D(view(B.data, :, 1), B.mesh[1], w)
            else
                newB.data[i, j] = CompositeGrids.Interp.linear2D(B.data, B.mesh[1], B.mesh[2], w, v)
            end
        end
    end
    return newA, newB
end

function symmetrize_b(B)
    wgrid = B.mesh[1]
    sB = similar(B)
    aB = similar(B)
    for i in 1:length(wgrid)
        for j in 1:length(wgrid)
            sB.data[i, j] = (B[i, j] + B[j, i]) / 2.0
            aB.data[i, j] = (B[i, j] - B[j, i]) / 2.0
        end
    end
    return sB, aB
end

function interp_AB(β, A, B, param;
    Nk=12, minterval=0.5π / β, order=6)
    # interp AB to different temperature
    α = 0.882
    Euv = A.mesh[1][end]
    wgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [α / β, Euv], [α / β, 0.1param.EF], Nk, minterval, order)
    newA = GreenFunc.MeshArray(wgrid; dtype=Float64)
    newB = GreenFunc.MeshArray(wgrid, wgrid; dtype=Float64)
    newA.data .= CompositeGrids.Interp.interp1DGrid(A.data, A.mesh[1], newA.mesh[1])
    for i in 1:length(wgrid)
        for j in 1:length(wgrid)
            newB.data[i, j] = CompositeGrids.Interp.linear2D(B.data, B.mesh[1], B.mesh[2], wgrid[i], wgrid[j])
        end
    end
    sB, aB = symmetrize_b(newB)
    newparam = ElectronGas.Parameter.Para(param; β=β)
    return newparam, newA, sB
end

function interp_AB_brutal(β, A, B, param; Ec=10param.EF)
    # interp AB to different temperature
    nec = floor(Int, Ec / 2 / π * β + 0.5)
    wgrid = GreenFunc.ImFreq(β, FERMION; grid=[i for i in 0:nec])
    oldwgrid = A.mesh[1]
    # println((oldwgrid[1], oldwgrid[end]))
    # println((wgrid[1], wgrid[end]))
    newA = GreenFunc.MeshArray(wgrid; dtype=Float64)
    newB = GreenFunc.MeshArray(wgrid, wgrid; dtype=Float64)
    newA.data .= CompositeGrids.Interp.interp1DGrid(A.data, A.mesh[1], [w for w in wgrid])
    for i in 1:length(wgrid)
        for j in 1:length(wgrid)
            newB.data[i, j] = CompositeGrids.Interp.linear2D(B.data, B.mesh[1], B.mesh[2], wgrid[i], wgrid[j])
        end
    end
    sB, aB = symmetrize_b(newB)
    newparam = ElectronGas.Parameter.Para(param; β=β)
    return newparam, newA, sB
end

function RS_inva(g, f, Q1, Q2)
    α = 0.882
    l = log(Q1 / α)
    L = log(Q2)
    return 1 + g * L + g * l * (1 - f * (1 + g * L))
end

function RS_interaction(w1, w2, param)
    # 2 from freq symmetry
    # 4π^2/kF^2*kF from momentum integral and Πs
    # π from definition of Vrs
    # 2/π from atan in Πs
    factor = 2 / param.kF * 4 * π^2 * π * 2 / π
    g, f, Ω, Ec = 0.2 * factor, 0.8, 0.2param.EF, 10param.EF
    if abs(w1) < Ω && abs(w2) < Ω
        return -g * (1 - f)
    elseif abs(w1) < Ec && abs(w2) < Ec
        return -g
    else
        return 0.0
    end
end

function RS_AB_brutal(β, A, B, param; Ec=10param.EF)
    # interp AB to different temperature
    nec = floor(Int, Ec / 2 / π * β + 0.5)
    wgrid = GreenFunc.ImFreq(β, FERMION; grid=[i for i in 0:nec])
    oldwgrid = A.mesh[1]
    newA = GreenFunc.MeshArray(wgrid; dtype=Float64)
    newB = GreenFunc.MeshArray(wgrid, wgrid; dtype=Float64)
    newA.data .= 1.0
    for i in 1:length(wgrid)
        for j in 1:length(wgrid)
            w, v = wgrid[i], wgrid[j]
            newB.data[i, j] = RS_interaction(w, v, param)
        end
    end
    sB, aB = symmetrize_b(newB)
    newparam = ElectronGas.Parameter.Para(param; β=β)
    return newparam, newA, sB
end

function Πs0(ωn, param; ω_c=0.1param.EF)
    @unpack me, β, μ, kF, EF = param
    return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
end

function Πs0wrapped(wgrid, param; ω_c=0.1param.EF)
    @unpack me, β, μ, kF, EF = param
    Π = GreenFunc.MeshArray(wgrid; dtype=Float64)
    Π.data .= [Πs0(w, param; ω_c=ω_c) for w in wgrid]
    return Π
end

function calcR!(R, A, B, Π, param)
    wgrid = A.mesh[1]
    result = similar(R.data)
    for (wi, w) in enumerate(wgrid)
        integrand = zeros(Float64, length(wgrid))
        for (vi, v) in enumerate(wgrid)
            integrand[vi] = B[wi, vi] * Π[vi] * R[vi]
        end
        # factor = 1 / 2 / π / (4 * π^2) * param.kF^2
        factor = 1 / 2 / π / (4 * π^2) * param.kF^2
        result[wi] = A[wi] + CompositeGrids.Interp.integrate1D(integrand, wgrid) * factor
    end
    R.data .= result
end

function calcR_brutal!(R, A, B, Π, param)
    wgrid = A.mesh[1]
    result = similar(R.data)
    result .= 0.0
    for (vi, v) in enumerate(wgrid)
        for (wi, w) in enumerate(wgrid)
            result[wi] += B[wi, vi] * Π[vi] * R[vi]
        end
    end
    # factor = param.kF^2 / (param.β * 4 * π^2)
    factor = param.kF^2 / (param.β * 4 * π^2)
    result .= result .* factor
    R.data .= result .+ A.data
end

function pcf_loop_ab(A, B, param; ω_c=0.1param.EF,
    α=0.9, Nmax=1e4)
    wgrid = A.mesh[1]
    Π = Πs0wrapped(wgrid, param; ω_c=ω_c)
    R = similar(A)
    R.data .= A.data
    Rsum = similar(R)
    Rsum.data .= R.data ./ (1 - α)
    iw0 = 1

    invR0 = 0.0
    diff = 1.0
    converge = false
    n = 1
    while (!converge && n < Nmax)
        calcR!(R, A, B, Π, param)
        Rsum.data .= Rsum.data .* α .+ R.data
        R.data .= Rsum.data .* (1 - α)
        converge = isapprox(1 / R[iw0], invR0, rtol=1e-10, atol=1e-10)
        invR0 = 1 / R[iw0]
        n = n + 1
        # println("invR0=$invR0")
    end

    return invR0, R

end

function pcf_loop_ab_brutal(A, B, param; ω_c=0.1param.EF,
    α=0.9, Nmax=1e4)
    wgrid = A.mesh[1]
    Π = Πs0wrapped(wgrid, param; ω_c=ω_c)
    R = similar(A)
    R.data .= A.data
    Rsum = similar(R)
    Rsum.data .= R.data ./ (1 - α)
    iw0 = 1

    invR0 = 0.0
    diff = 1.0
    converge = false
    n = 1
    while (!converge && n < Nmax)
        calcR_brutal!(R, A, B, Π, param)
        Rsum.data .= Rsum.data .* α .+ R.data
        R.data .= Rsum.data .* (1 - α)
        converge = isapprox(1 / R[iw0], invR0, rtol=1e-10, atol=1e-10)
        invR0 = 1 / R[iw0]
        n = n + 1
        # println("invR0=$invR0")
        # println((R[1], R[end-1]))
    end

    return invR0, R

end

using Test
@testset "pcf resummed" begin
    # fname = "run/data/PCFresumdlr_3000010.jld2"
    fname = "run/data/PCFresumrs3_3000022.jld2"
    param, A, B = load_AB(fname)
    # println(size(A))
    # println(size(B))
    println((B[1, 1], B[1, end]))

    fname = "run/data/PCFresumdlr_3000044.jld2"
    param, B = load_B(fname)
    println((B[1, 1], B[1, end]))

    # println(size(B))
    # println(param)
    # println((A[1], A[end]))
    # println((B[1, 1], B[1, end], B[end, 1], B[end, end]))
    A, B = extend_AB(A, B, param)
    num = 5
    betas = [200 * 2^(i - 1) for i in 1:num]
    lamus = zeros(Float64, length(betas))
    for i in 1:length(betas)
        beta = betas[i]
        # newparam, newA, newB = interp_AB(beta / param.EF, A, B, param)
        newparam, newA, newB = interp_AB_brutal(beta / param.EF, A, B, param)
        # newparam, newA, newB = RS_AB_brutal(beta / param.EF, A, B, param)

        newB.data .*= newparam.kF
        # println(newparam.β)
        # println((newA[1], newA[end]))
        # println((newB[1, 1], newB[1, end], newB[end, 1], newB[end, end]))

        # lamu, R = pcf_loop_ab(newA, newB, newparam)
        # lamu, R = pcf_loop_ab_brutal(newA, newB, newparam; ω_c=40param.EF)
        lamu, R = pcf_loop_ab_brutal(newA, newB, newparam)
        println("lamu=$lamu")
        lamus[i] = lamu
    end
end
