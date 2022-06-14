"""
calculate the Z factor of UEG using MC
"""

using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random
using MCIntegration
using Lehmann

const steps = 1e7
# include("parameter.jl")
dim = 3
beta = 100.0
rs = 5.0
mass2 = 1.e-6
# const para = Parameter.rydbergUnit(1.0 / beta, rs, dim, Λs = mass2)
const para = Parameter.defaultUnit(1.0 / beta, rs, dim, Λs=mass2)
const kF = para.kF
const EF = para.EF
const β = para.β

# qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 *kF , 8)
qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 1e-6 * kF, 8)
τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

# dlr = DLRGrid(Euv = 10EF, β = β, rtol = 1e-10, isFermi = false, symmetry = :ph)
# W = zeros(length(qgrid), dlr.size)
# for (qi, q) in enumerate(qgrid.grid)
#     for (ni, n) in enumerate(dlr.n)
#         # W[qi, ni] = Interaction.RPA(q, n, para, regular = true, pifunc = Polarization.Polarization0_FiniteTemp)[1]
#         W[qi, ni] = Interaction.RPA(q, n, para, regular = true, pifunc = Polarization.Polarization0_ZeroTemp)[1]
#     end
# end
# const dW0 = matfreq2tau(dlr, W, τgrid.grid, axis = 2)
Ws, Wa = Interaction.RPAwrapped(10EF, 1e-10, qgrid.grid, para)
const dW0 = real.(matfreq2tau(Ws.dlrGrid, Ws.dynamic, τgrid.grid, axis=4))[1, 1, :, :]

function integrand(config)
    if config.curr == 1
        return eval2(config)
    else
        error("impossible!")
    end
end

function interactionDynamic(qd, τIn, τOut)
    dτ = abs(τOut - τIn)
    kDiQ = sqrt(dot(qd, qd))

    # vd = 4π * para.e0^2 / (kDiQ^2 + para.Λs)
    vs, va = Interaction.coulomb_2d(kDiQ, para)
    # vs, va = Interaction.coulomb(kDiQ, para)
    if kDiQ <= qgrid.grid[1]
        q = qgrid.grid[1] + 1.0e-6
        wd = vs * Interp.linear2D(dW0, qgrid, τgrid, q, dτ)
        # the current interpolation vanishes at q=0, which needs to be corrected!
    else
        wd = vs * Interp.linear2D(dW0, qgrid, τgrid, kDiQ, dτ) # dynamic interaction, don't forget the singular factor vq
    end

    return wd
end

function eval2(config)

    K, T = config.var[1], config.var[2]
    k, τ = K[1], T[1]
    k0 = zeros(para.dim)
    k0[end] = kF # external momentum
    kq = k + k0
    ω = (dot(kq, kq) - kF^2) / (2 * para.me)
    g = Spectral.kernelFermiT(τ, ω, β)
    dW = interactionDynamic(k, 0.0, τ)
    phase = 1.0 / (2π)^para.dim
    return g * dW * phase
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    τ = config.var[2][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[1] += weight * sin(π / β * τ) / abs(weight) * factor
        config.observable[2] += weight * sin(3π / β * τ) / abs(weight) * factor
    else
        return
    end
end

function fock(extn)
    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Tau(β, β / 2.0)

    dof = [[1, 1],] # degrees of freedom of the Fock diagram
    obs = zeros(2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T), dof, obs; para=para)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=16)

    if isnothing(avg) == false
        @printf("%10.6f   %10.6f ± %10.6f\n", -1.0, avg[1], std[1])
        @printf("%10.6f   %10.6f ± %10.6f\n", 0.0, avg[2], std[2])

        dS_dw = (avg[1] - avg[2]) / (2π / β)
        error = (std[1] + std[2]) / (2π / β)
        println("dΣ/diω= $dS_dw ± $error")
        Z = (1 / (1 + dS_dw))
        Zerror = error / Z^2
        println("Z=  $Z ± $Zerror")
        # TODO: add errorbar estimation
        # println
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # using Gaston
    ngrid = [-1, 0, 1]
    @time fock(ngrid)
end