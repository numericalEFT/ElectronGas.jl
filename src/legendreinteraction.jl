"""
Provide Legendre decomposed interaction that could be useful for calculating self-energy and solving gap-function
equation of superconductivity problems.
"""
module LegendreInteraction

using Parameters, GreenFunc, Lehmann, LegendrePolynomials, CompositeGrids

export DCKernel

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

include(srcdir*"/parameter.jl")
using .Parameter

include(srcdir*"/convention.jl")
using .Convention

include(srcdir*"/polarization.jl")
using .Polarization

include(srcdir*"/interaction.jl")
using .Interaction

function interaction_dynamic(q, n, param, int_type, spin_state)
    # a wrapper for dynamic part of effective interaction
    # for rpa simply return rpa
    # for ko return ks+ka for singlet, ks-3ka for triplet
    if int_type == :rpa
        return RPA(q, n, param)
    elseif int_type == :ko
        ks, ka = KO(q, n, param)
        if spin_state==:singlet
            spin_factor = 1.0
        elseif spin_state==:triplet
            spin_factor = -3.0
        elseif spin_state==:sigma
            spin_factor = 3.0
        else
            throw(UndefVarError(spin_state))
        end
        return ks + spin_factor * ka
    else
        throw(UndefVarError(int_type))
    end
end

@inline function kernel_integrand(k, p, q, n, channel, param, int_type, spin_state)
    legendre_x = (k^2 + p^2 - q^2)/2/k/p
    if(abs(abs(legendre_x)-1)<1e-12)
        legendre_x = sign(legendre_x)*1
    end
    return q*Pl(legendre_x, channel)*V_Bare(q, param)*interaction_dynamic(q, n, param, int_type, spin_state)
end

@inline function kernel0_integrand(k, p, q, channel, param)
    legendre_x = (k^2 + p^2 - q^2)/2/k/p
    if(abs(abs(legendre_x)-1)<1e-12)
        legendre_x = sign(legendre_x)*1
    end
    @assert -1<=legendre_x<=1 "k=$k,p=$p,q=$q"
    return q*Pl(legendre_x, channel)*V_Bare(q, param)
end


struct DCKernel
    int_type::Symbol
    spin_state::Symbol
    channel::Int

    param::Parameter.Para

    kgrid::CompositeGrid.Composite
    qgrids::Vector{CompositeGrid.Composite}
    dlrGrid::DLRGrid

    kernel_bare::Array{Float64,2}
    kernel::Array{Float64,3}

    function DCKernel(param, Euv, rtol, Nk, maxK, minK, order, int_type, channel, spin_state=:auto)
        @unpack kF, beta = param
        EPS = 1e-16

        if spin_state==:sigma
            # for self-energy, always use ℓ=0
            channel = 0
        elseif spin_state==:auto
            # automatically assign spin_state, triplet for even, singlet for odd channel
            spin_state = (channel%2==0) ? (:triplet) : (:singlet)
        end

        bdlr = DLRGrid(Euv, beta, rtol, false, :ph)
        kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order )
        #println(kgrid.grid)
        qgrids = [CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [k, kF], Nk, minK, order) for k in kgrid.grid]
        qgridmax = maximum([qg.size for qg in qgrids])
        #println(qgridmax)

        kernel_bare = zeros(Float64, (length(kgrid.grid), (qgridmax)))
        kernel = zeros(Float64, (length(kgrid.grid), (qgridmax), length(bdlr.n)))

        int_grid_base = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.1*maxK], [0.0, 2kF], 2Nk, 0.01minK, 2)
        for (ki, k) in enumerate(kgrid.grid)
            for (pi, p) in enumerate(qgrids[ki].grid)
                if abs(k - p) > EPS

                    kmp = abs(k-p)<EPS ? EPS : abs(k-p)
                    kpp = k + p
                    im, ip = floor(int_grid_base, kmp), floor(int_grid_base, kpp)
                    int_panel = Float64[]

                    push!(int_panel, kmp)
                    if im<ip
                        for i in im+1:ip
                            push!(int_panel, int_grid_base[i])
                        end
                    end
                    push!(int_panel, kpp)

                    int_panel = SimpleGrid.Arbitrary{Float64}(int_panel)
                    SubGridType = SimpleGrid.GaussLegendre{Float64}
                    subgrids = subgrids = Vector{SubGridType}([])
                    for i in 1:int_panel.size-1
                        _bound = [int_panel[i],int_panel[i+1]]
                        push!(subgrids, SubGridType(_bound,order))
                    end
                    
                    int_grid=CompositeGrid.Composite{Float64,SimpleGrid.Arbitrary{Float64},SubGridType}(int_panel,subgrids)

                    data = [kernel0_integrand(k, p, q, channel, param) for q in int_grid.grid]
                    kernel_bare[ki, pi] = Interp.integrate1D(data, int_grid)

                    for (ni, n) in enumerate(bdlr.n)
                        data = [kernel_integrand(k,p,q,n,channel,param,int_type,spin_state) for q in int_grid.grid]
                        kernel[ki, pi, ni] = Interp.integrate1D(data, int_grid)
                        @assert isfinite(kernel[ki,pi,ni]) "fail kernel at $ki,$pi,$ni, with $(kernel[ki,pi,ni])"
                    end

                else
                    kernel_bare[ki,pi] = 0
                    for (ni, n) in enumerate(bdlr.n)
                        kernel[ki,pi,ni] = 0
                    end
                end
            end
        end
        
        return new(int_type, spin_state, channel, param, kgrid, qgrids, bdlr, kernel_bare, kernel)
    end

    # function DCKernel(fineKernel::DCKernel, beta)
    #     # extrapolate to kernel of β from fineKernel of lower temperature
    #     @assert beta<fineKernel.param.beta "can only extrapolate from low temp to high temp!"
    #     int_type, spin_state, channel = fineKernel.int_type, fineKernel.spin_state, fineKernel.channel
    #     param = Parameter.Para(fineKernel.param, beta=beta) # new param with new beta

    #     Euv, rtol = fineKernel.dlr.Euv, fineKernel.dlr.rtol
    #     bdlr = DLRGrid(Euv, beta, rtol, false, :ph)

    #     kgrid, qgrids = fineKernel.kgrid, fineKernel.qgrids

    #     kernel_bare = fineKernel.kernel_bare

    #     fineKernel_dlr = real(matfreq2dlr(fineKernel.dlr, fineKernel.kernel; axis=3))
    #     kernel = real(dlr2matfreq(fineKernel.dlr, fineKernel_dlr, bdlr.n ./ (beta/fineKernel.param.beta); axis=3))

    #     return new(int_type, spin_state, channel, param, kgrid, qgrids, bdlr, kernel_bare, kernel)
    # end

end

end

if abspath(PROGRAM_FILE) == @__FILE__
    # using Plots

    param = LegendreInteraction.Parameter.defaultUnit(1000.0, 1.0)

    kernel = LegendreInteraction.DCKernel(param, 100*param.EF, 1e-8, 5, 10*param.kF, 1e-7*param.kF, 4, :rpa,0,:sigma)
    # kernel2 = LegendreInteraction.DCKernel(kernel, 500.0)
    kF_label = searchsortedfirst(kernel.kgrid.grid, kernel.param.kF)
    qF_label = searchsortedfirst(kernel.qgrids[kF_label].grid, kernel.param.kF)
    println(kernel.kernel[kF_label,qF_label,:])
    # println(kernel2.kernel[kF_label,qF_label,:])

    # p = plot(kernel.dlr.ωn[1:8], kernel.kernel[kF_label,qF_label,1:8])
    # plot!(p, kernel2.dlr.ωn[1:8], kernel2.kernel[kF_label,qF_label,1:8])
    # display(p)
    # readline()

end
