using ElectronGas
using JLD2
using ElectronGas.CompositeGrids
using ElectronGas.Parameters
using ElectronGas.GreenFunc
using LsqFit
using DelimitedFiles

function load_zph(fname)
    f = jldopen(fname, "r")
    n, Z = f["n"], f["Z"]
    return SimpleG.Arbitrary(n), Z
end

# ngrid, Zph = load_zph("run/data/sigma.jld2")
# ngrid, Zph = load_zph("run/data/sigma_12800.0.jld2")
# ngrid, Zph = load_zph("run/data/sigma_51200_0005.jld2")
# ngrid, Zph = load_zph("run/data/sigma_51200_00025.jld2")
# ngrid, Zph = load_zph("run/data/sigma_51200_00025_nokf.jld2")
# ngrid, Zph = load_zph("run/data/sigma_51200_0005_nokf.jld2")
ngrid, Zph = load_zph("run/data/sigma_51200_191916_0005.jld2")
# ngrid, Zph = nothing, nothing
# Zph .= 1.0
# println(ngrid, Zph)
betaZ = 51200

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
    Nk=24, minterval=0.5π / β, order=6)
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

function interp_AB_brutal_step(β, A, B, param; Ec=0.1param.EF)
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

# function Πs0(ωn, param; ω_c=0.1param.EF, ngrid=nothing, Z=nothing)
function Πs0(ωn, param; ω_c=0.1param.EF, ngrid=ngrid, Zph=Zph)
    @unpack me, β, μ, kF, EF = param
    # return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
    if isnothing(ngrid)
        if abs(ωn) > ω_c
            return 0.0
        else
            # return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
            return π * me / (kF) / abs(ωn)
        end
    else
        if abs(ωn) > ω_c
            return 0.0
        else
            n = (ωn / π * β / param.beta * betaZ - 1) / 2
            zph = Interp.interp1D(Zph, ngrid, n)
            ωn = ωn / zph
            return π * me / (kF) / abs(ωn)
            # return 2.0 * me / (kF) / abs(ωn) * atan(ω_c / abs(ωn))
        end
    end
end

function Πs0wrapped(wgrid, param; ω_c=0.1param.EF)
    @unpack me, β, μ, kF, EF = param
    Π = GreenFunc.MeshArray(wgrid; dtype=Float64)
    Π.data .= [Πs0(w, param; ω_c=ω_c) for w in wgrid]
    return Π
end

function Πstail(wi, wf, param; ω_c=0.1param.EF)
    result = 0.0
    ni = floor(Int, wi / 2 / π * param.β + 0.5)
    nf = floor(Int, wf / 2 / π * param.β + 0.5)
    for i in ni:nf
        ω = π / param.β * (2 * i + 1)
        result += Πs0(ω, param; ω_c=ω_c)
    end
    return result
end

function calcR!(R, A, B, Π, param; tail=0.0)
    wgrid = A.mesh[1]
    result = similar(R.data)
    for (wi, w) in enumerate(wgrid)
        integrand = zeros(Float64, length(wgrid))
        for (vi, v) in enumerate(wgrid)
            integrand[vi] = B[wi, vi] * Π[vi] * R[vi]
        end
        # factor = 1 / 2 / π / (4 * π^2) * param.kF^2
        integraltail = B[wi, end] * tail * R[end]
        factor = 1 / 2 / π / (4 * π^2) * param.kF^2
        tailfactor = param.kF^2 / (param.β * 4 * π^2)
        result[wi] = A[wi] + CompositeGrids.Interp.integrate1D(integrand, wgrid) * factor + integraltail * tailfactor
    end
    R.data .= result
end

function calcR_brutal!(R, A, B, Π, param; tail=0.0)
    wgrid = A.mesh[1]
    result = similar(R.data)
    result .= 0.0
    for (wi, w) in enumerate(wgrid)
        for (vi, v) in enumerate(wgrid)
            result[wi] += B[wi, vi] * Π[vi] * R[vi]
        end
        result[wi] += B[wi, end] * tail * R[end]
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
    tail = Πstail(wgrid[end], 10.0 * wgrid[end], param; ω_c=ω_c)
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
        calcR!(R, A, B, Π, param; tail=tail)
        Rsum.data .= Rsum.data .* α .+ R.data
        R.data .= Rsum.data .* (1 - α)
        converge = isapprox(1 / R[iw0], invR0, rtol=1e-10, atol=1e-10)
        invR0 = 1 / R[iw0]
        n = n + 1
        # println("invR0=$invR0")
    end

    return invR0, R

end

function pcf_loop_ab_brutal_step(A, B, param; ω_c=0.1param.EF,
    α=0.9, Nmax=1e4)
    wgrid = A.mesh[1]
    Π = Πs0wrapped(wgrid, param; ω_c=ω_c)
    # tail = Πstail(wgrid[end], 200.0 * wgrid[end], param; ω_c=ω_c)
    tail = 0.0
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
        calcR_brutal!(R, A, B, Π, param; tail=tail)
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

@. model(x, p) = p[1] * x + p[2]
function fit_invR0(invR0, lnbetas)
    p0 = [1.0, 1.0]
    fit = curve_fit(model, lnbetas, invR0, p0)
    return fit
end

function crit_beta(betas, lamus; init=0, fin=length(betas))
    betas = betas[init:fin]
    lamus = lamus[init:fin]
    lnbetas = log10.(betas)
    fitresult = fit_invR0(lamus, lnbetas)
    fitp = coef(fitresult)
    println(fitp)
    return -fitp[2] / fitp[1]
end

function avg_phonon(n, param)
    # compute phonon averaged on fermi surface
    xgrid = SimpleG.GaussLegendre{Float64}([-1.0, 1.0], 20)
    data = zeros(Float64, length(xgrid))
    for i in 1:length(xgrid)
        data[i] = Interaction.phonon(n, param.kF * sqrt(2 - 2 * xgrid[i]), param)[1]
    end
    return Interp.integrate1D(data, xgrid)
    # return Interp.integrate1D(data, xgrid) / π^2 * 4
end

function add_phonon(B, param)
    # add phonon to B
    wgrid = B.mesh[1]
    for (wi, w) in enumerate(wgrid)
        for (vi, v) in enumerate(wgrid)
            dw = abs(w - v)
            n = dw / 2 / π * param.β
            B.data[wi, vi] -= avg_phonon(n, param)
            dw = abs(w + v)
            n = dw / 2 / π * param.β
            B.data[wi, vi] -= avg_phonon(n, param)
        end
    end
    return B
end

# fname = "run/data/PCFresumdlr_3000800001.jld2"
# f = jldopen(fname, "r")
# param, B0 = f["param"], f["B"]
# println(param)

# param = ElectronGas.Parameter.rydbergUnit(1 / 23456, 1.0, 3)
param = ElectronGas.Parameter.rydbergUnit(1 / 23456, 1.91916, 3)
wgrid = CompositeGrids.CompositeG.LogDensedGrid(:cheb, [0.882 / param.β, 10param.EF], [0.882 / param.β, 0.1param.EF], 64, 0.882 / param.β, 10)
B0 = GreenFunc.MeshArray(wgrid, wgrid; dtype=Float64)
# B0.data .= -4.41997
# B0.data .= -4.135597143894588
# B0.data .= 0.6770001354095592
B0.data .= 0.0
# println(B0[1] * param.kF^2)
# println(param.kF^2)
println(param)
wgrid = B0.mesh[1]
A = MeshArray(wgrid; dtype=Float64)
B = MeshArray(wgrid, wgrid; dtype=Float64)
A.data .= 1
B.data .= B0[1]

param = ElectronGas.Parameter.Para(param; eph=0.4 * 4 * π^2, ω_D=0.005param.EF)
B = add_phonon(B, param)

A, B = extend_AB(A, B, param)
# num = 7
# betas = [400 * 2^(i - 1) for i in 1:num]
# num = 24
# betas = [50 * sqrt(2)^(i - 1) for i in 1:num]
num = 37
betas = [50 * 2^((i - 1) / 4) for i in 1:num]
lamus = zeros(Float64, length(betas))
for i in 1:length(betas)
    beta = betas[i]
    newparam, newA, newB = interp_AB(beta / param.EF, A, B, param)
    # newparam, newA, newB = interp_AB_brutal_step(beta / param.EF, A, B, param)
    # newA.data .= 1.0
    # newB.data .= B0[1]
    # newparam, newA, newB = RS_AB_brutal(beta / param.EF, A, B, param)

    # println(newparam.β)
    # println((newA[1], newA[end]))
    # println((newB[1, 1], newB[1, end], newB[end, 1], newB[end, end]))

    lamu, R = pcf_loop_ab(newA, newB, newparam)
    # lamu, R = pcf_loop_ab_brutal(newA, newB, newparam; ω_c=40param.EF)
    # lamu, R = pcf_loop_ab_brutal_step(newA, newB, newparam)
    println("β=$beta, lamu=$lamu")
    lamus[i] = lamu
end
log10tc = -crit_beta(betas, lamus; init=16)
# println("$log10tc, Tc=$(10^log10tc)")
# log10tc = -crit_beta(betas, lamus; init=10, fin=7)
println("$log10tc, Tc=$(10^log10tc)")

# savefname = "./run/rpcf3D_phrpaPiph2_rs1.0_l0_vlarge0.txt"
# savefname = "./run/rpcf3D_phrpachi_rs1.0_l0_vlarge0.txt"
# savefname = "./run/rpcf3D_phrpaPiph_rs1.0_l0_vlarge0.txt"

# savefname = "./run/rpcf3D_rpaPiph_rs1.91916_l0_vlarge0.txt"
# savefname = "./run/rpcf3D_ph1rpa_rs1.91916_l0_vlarge0.txt"
# savefname = "./run/rpcf3D_ph1_rs1.91916_l0_vlarge0.txt"
savefname = "./run/rpcf3D_ph1_rs1.91916_l0_sigma.txt"

open(savefname, "w") do io
end
for i in 1:num
    data = [betas[i] lamus[i] lamus[i] 0 1.0]
    open(savefname, "a+") do io
        writedlm(io, data, ' ')
    end
end
