"""
Calculate gap-function equation
"""

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.Parameters
using ElectronGas.Lehmann
using ElectronGas.CompositeGrids

function G02wrapped(fdlr, kgrid, param)
    # return G(K)G(-K)
    @unpack me, kF, β, EF = param

    green_dyn = zeros(Float64, (kgrid.size, fdlr.size))
    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(fdlr.n)
            ω = k^2/2/me
            green_dyn[ki,ni] = 1/(
                ( (2n + 1) * π / β) ^ 2
                + (ω) ^ 2
            )
        end
    end

    return green_dyn
end

function Δinit(fdlr, kgrid)
    @unpack me, kF, β, EF = param

    delta = zeros(Float64, (kgrid.size, fdlr.size))
    delta0 = zeros(Float64, kgrid.size)

    for (ki, k) in enumerate(kgrid)
        delta0[ki] = 1.0
        for (τi, τ) in enumerate(fdlr.τ)
            delta[ki, τi] = 1.0
        end
    end

    return delta, delta0
end

function normalizeΔ()
end

function calcF!()
end

function calcΔ!()
end

function gapIteration(param, fdlr, kgrid, qgrids,  kernel, kernel_bare, G2, Δ, Δ0;
                      Nstep = 1e3, rtol = 1e-3, shift = 2.0)

end



if abspath(PROGRAM_FILE) == @__FILE__

    #--- parameters ---


    param = Parameter.defaultUnit(1/1000.0, 3.0)
    Euv, rtol = 100*param.EF, 1e-10
    maxK, minK = 10param.kF, 1e-6param.kF
    Nk, order = 8, 4
    int_type = :rpa

    #--- prepare kernel ---
    W = LegendreInteraction.DCKernel0(param;
                                      Euv = Euv, rtol = rtol, Nk = Nk, maxK = maxK, minK = minK, order = order,
                                      int_type = int_type)

    fdlr = Lehmann.DLRGrid(Euv, param.β, rtol, true, :ph)
    bdlr = W.dlrGrid
    kgrid = W.kgrid
    qgrids = W.qgrids

    kernel_bare = W.kernel_bare
    kernel_freq = W.kernel
    kernel = real(Lehmann.matfreq2tau(bdlr, kernel_freq, fdlr.τ, bdlr.n; axis = 3))

    kF_label = searchsortedfirst(kgrid.grid, param.kF)
    qF_label = searchsortedfirst(qgrids[kF_label].grid, param.kF)

    println("dynamic kernel at (kF, kF):")
    print(view(kernel, kF_label, qF_label, :))

    #--- prepare G2 ---

    G2 = G02wrapped(fdlr, kgrid, param)

    println("G2 at kF:")
    print(view(G2, kF_label, :))

    #--- prepare Δ ---
    Δ, Δ0 = Δinit(fdlr, kgrid)


    # Δ = gapIteration(param, 10, Euv, rtol, Nk, 10*param.kF, 1e-7*param.kF, order, :rpa)

    # kgrid = Δ.spaceGrid
    # kF = kgrid.panel[3]
    # kF_label = searchsortedfirst(kgrid.grid, kF)

    # println(Δ.instant[1,1,:])
    # println(Δ.dynamic[1,1,kF_label,:])
end
