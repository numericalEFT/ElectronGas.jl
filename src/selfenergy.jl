"""
Calculate self-energy
"""

module SelfEnergy

using ..Parameter, ..Convention, ..Polarization, ..Interaction, ..LegendreInteraction
using ..Parameters, ..GreenFunc, ..Lehmann, ..LegendrePolynomials, ..CompositeGrids

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

function G0wrapped(Euv,rtol,sgrid,param)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0, EF = param

    green = GreenFunc.Green2DLR{ComplexF64}(:g0,GreenFunc.IMFREQ,beta,true,Euv,sgrid,1)
    green_dyn = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    green_ins = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1,1,ki,ni] = 1/(im*(π/beta*(2n+1)) - (k^2/2/me-EF) )
        end
    end
    green.dynamic=green_dyn
    return green
end

# function calcΣ(kernal, kernal_bare, fdlr, kgrid, qgrids)
function calcΣ(G::GreenFunc.Green2DLR, W::LegendreInteraction.DCKernel)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0 = W.param

    kgrid = W.kgrid
    qgrids = W.qgrids
    fdlr = G.dlrGrid
    bdlr = W.dlrGrid
    G=GreenFunc.toTau(G)

    # prepare kernel, interpolate into τ-space with fdlr.τ
    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis = 3)

    # container of Σ
    Σ = GreenFunc.Green2DLR{ComplexF64}(:sigma, GreenFunc.IMTIME,beta,true,fdlr.Euv,kgrid,1)
    Σ_ins = zeros(ComplexF64, (1,1,length(kgrid.grid)))
    Σ_dyn = zeros(ComplexF64, (1,1,length(kgrid.grid), fdlr.size))

    # tbc
    G_ins = tau2tau(G.dlrGrid, G.dynamic, [beta, ], G.timeGrid.grid; axis = 4)[1,1,:,1]

    for (ki, k) in enumerate(kgrid.grid)

        for (τi, τ) in enumerate(fdlr.τ)
            Gq = CompositeGrids.Interp.interp1DGrid(G.dynamic[1,1,:,τi],kgrid,qgrids[ki].grid)
            integrand = kernel[ki, 1:qgrids[ki].size, τi] .* Gq ./ k .* qgrids[ki].grid
            Σ_dyn[1,1,ki, τi] = CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
            @assert isfinite(Σ_dyn[1,1,ki, τi]) "fail Δ at $ki, $τi"

            if τi == 1
                Gq = CompositeGrids.Interp.interp1DGrid(G_ins,kgrid,qgrids[ki].grid)
                integrand = kernel_bare[ki, 1:qgrids[ki].size] .* Gq ./ k .* qgrids[ki].grid
                Σ_ins[1,1,ki] += CompositeGrids.Interp.integrate1D(integrand, qgrids[ki])
                @assert isfinite(Σ_ins[1,1,ki]) "fail Δ0 at $ki"
            end

        end
    end

    Σ.dynamic, Σ.instant = Σ_dyn./(4*π^2), Σ_ins./(4*π^2)
    return Σ
end

end

if abspath(PROGRAM_FILE) == @__FILE__

    param = SelfEnergy.LegendreInteraction.Parameter.defaultUnit(1000.0, 1.0)
    Euv, rtol = 100*param.EF, 1e-8
    kernel = SelfEnergy.LegendreInteraction.DCKernel(param, Euv, rtol, 5, 10*param.kF, 1e-7*param.kF, 4, :rpa,0,:sigma)

    G0 = SelfEnergy.G0wrapped(Euv, rtol, kernel.kgrid, param)

    Σ = SelfEnergy.calcΣ(G0, kernel)

    println(Σ.dynamic)

end