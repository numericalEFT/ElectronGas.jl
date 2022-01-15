"""
Calculate self-energy
"""

module SelfEnergy

using ..Parameter, ..Convention, ..Polarization, ..Interaction, ..LegendreInteraction
using ..Parameters, ..GreenFunc, ..Lehmann, ..LegendrePolynomials, ..CompositeGrids

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

function G0wrapped(Euv,rtol,sgrid,param)
    @unpack me, kF, beta, EF = param

    green = GreenFunc.Green2DLR{ComplexF64}(:g0,GreenFunc.IMFREQ,beta,true,Euv,sgrid,1)
    green_dyn = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1,1,ki,ni] = 1/(im*(π/beta*(2n+1)) - (k^2/2/me-EF) )
        end
    end
    green.dynamic=green_dyn
    return green
end

function Gwrapped(Σ::GreenFunc.Green2DLR, param)
    @unpack me, kF, beta, EF = param
    Σ_freq = GreenFunc.toMatFreq(Σ)
    green =  Green2DLR{ComplexF64}(
        :G, GreenFunc.IMFREQ,Σ_freq.β, Σ_freq.isFermi, Σ_freq.dlrGrid.Euv, Σ_freq.spaceGrid, Σ_freq.color;
        timeSymmetry = Σ_freq.timeSymmetry, rtol = Σ_freq.dlrGrid.rtol)

    green_dyn = zeros(ComplexF64, (green.color, green.color, green.spaceGrid.size, green.timeGrid.size))
    for (ki, k) in enumerate(green.spaceGrid)
        for (ni, n) in enumerate(green.dlrGrid.n)
            green_dyn[1,1,ki,ni] = 1/(im*(π/beta*(2n+1)) - (k^2/2/me-EF) + Σ.dynamic[1,1,ki,ni] + Σ.instant[1,1,ki])
        end
    end
    green.dynamic=green_dyn
    return green
end

# function calcΣ(kernal, kernal_bare, fdlr, kgrid, qgrids)
function calcΣ(G::GreenFunc.Green2DLR, W::LegendreInteraction.DCKernel)
    @unpack beta= W.param

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

function G0W0(param, Euv, rtol, Nk, maxK, minK, order, int_type)
    kernel = SelfEnergy.LegendreInteraction.DCKernel0(param;
                                                      Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order, int_type=int_type, spin_state=:sigma)
    G0 = G0wrapped(Euv, rtol, kernel.kgrid, param)
    Σ = calcΣ(G0, kernel)

    return Σ
end

function zfactor(Σ::GreenFunc.Green2DLR)
    kgrid = Σ.spaceGrid
    kF = kgrid.panel[3]
    beta = Σ.dlrGrid.β

    println("kF=$kF")
    kF_label = searchsortedfirst(kgrid.grid, kF)
    Σ_freq = GreenFunc.toMatFreq(Σ, [0,1])

    ΣI = imag(Σ_freq.dynamic[1,1,kF_label,:])

    # for correct sign of ΣI should be 1/(1 - (ΣI[2]-ΣI[1])/2/π*beta)
    Z0 = 1/(1 + (ΣI[2]-ΣI[1])/2/π*beta)
    return Z0
end

end
