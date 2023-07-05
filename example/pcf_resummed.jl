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
    Nk=8, minterval=π / β, order=4)
    # interp AB to different temperature
    α = 0.882
    Euv = A.mesh[1][end]
    wgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [α / β, Euv], [α / β, 0.1param.EF], Nk, minterval, order)
    oldwgrid = A.mesh[1]
    newA = GreenFunc.MeshArray(wgrid; dtype=Float64)
    newB = GreenFunc.MeshArray(wgrid, wgrid; dtype=Float64)
    newA.data .= CompositeGrids.Interp.interp1DGrid(A.data, A.mesh[1], newA.mesh[1])
    for i in 1:length(wgrid)
        for j in 1:length(wgrid)
            newB.data[i, j] = CompositeGrids.Interp.linear2D(B.data, B.mesh[1], B.mesh[2], wgrid[i], wgrid[j])
        end
    end
    sB, aB = symmetrize_b(newB)
    return newA, sB
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

function calcR!(R, A, B, Π)
    wgrid = A.mesh[1]
    result = similar(R.data)
    for (wi, w) in enumerate(wgrid)
        integrand = zeros(Float64, length(wgrid))
        for (vi, v) in enumerate(wgrid)
            integrand[vi] = B[wi, vi] * Π[vi] * R[vi]
        end
        factor = 1 / 2π / (4π^2)
        result[wi] = A[wi] + CompositeGrids.Interp.integrate1D(integrand, wgrid) * factor
    end
    R.data .= result
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
        calcR!(R, A, B, Π)
        Rsum.data .= Rsum.data .* α .+ R.data
        R.data .= Rsum.data .* (1 - α)
        converge = isapprox(1 / R[iw0], invR0, rtol=1e-10, atol=1e-10)
        invR0 = 1 / R[iw0]
        n = n + 1
        # println("invR0=$invR0")
    end

    return invR0, R

end

using Test
@testset "pcf resummed" begin
    fname = "run/data/PCFresumdlr_3000010.jld2"
    param, A, B = load_AB(fname)
    A, B = extend_AB(A, B, param)
    num = 11
    betas = [400 * 2^(i - 1) for i in 1:num]
    lamus = zeros(Float64, length(betas))
    for i in 1:length(betas)
        beta = betas[i]
        newA, newB = interp_AB(beta / param.EF, A, B, param)
        lamu, R = pcf_loop_ab(newA, newB, param)
        println("lamu=$lamu")
        lamus[i] = lamu
    end
end
