"""
Provide Legendre decomposed interaction that could be useful for calculating self-energy and solving gap-function
equation of superconductivity problems.
"""
module LegendreInteraction

# using Parameters, GreenFunc, Lehmann, LegendrePolynomials, CompositeGrids
# include(srcdir*"/parameter.jl")
# using .Parameter
# include(srcdir*"/convention.jl")
# using .Convention
# include(srcdir*"/polarization.jl")
# using .Polarization
# include(srcdir*"/interaction.jl")
# using .Interaction

using ..Parameter, ..Convention, ..Polarization, ..Interaction
using ..Parameters, ..GreenFunc, ..Lehmann, ..LegendrePolynomials, ..CompositeGrids

export DCKernel

@inline function spin_factor(spin_state)
    if spin_state == :singlet
        factor = 1.0
    elseif spin_state == :triplet
        factor = -3.0
    elseif spin_state == :sigma
        factor = 3.0
    else
        throw(UndefVarError(spin_state))
    end
    return factor
end

function interaction_dynamic(q, n, param, int_type, spin_state)
    # a wrapper for dynamic part of effective interaction
    # for rpa simply return rpa
    # for ko return ks+ka for singlet, ks-3ka for triplet
    @unpack dim = param
    if dim != 2 && dim != 3
        throw(UndefVarError(dim))
    end

    if int_type == :rpa
        if dim == 3
            ks, ka = RPA(q, n, param; regular = true)
        elseif dim == 2
            ks, ka = RPA(q, n, param; Vinv_Bare = Interaction.coulombinv_2d)
        end
    elseif int_type == :ko
        if dim == 3
            ks, ka = KO(q, n, param)
        elseif dim == 2
            ks, ka = KO(q, n, param; Vinv_Bare = Interaction.coulombinv_2d)
        end
    else
        throw(UndefVarError(int_type))
    end

    return ks + spin_factor(spin_state) * ka
end

@inline function interaction_instant(q, param, spin_state)
    @unpack dim = param
    if dim != 2 && dim != 3
        throw(UndefVarError(dim))
    end

    if dim == 3
        Vs, Va = coulomb(q, param)
    elseif dim == 2
        Vs, Va = coulomb_2d(q, param)
    end

    return (Vs + spin_factor(spin_state) * Va)
end

@inline function kernel_integrand(k, p, q, n, channel, param, int_type, spin_state)
    legendre_x = (k^2 + p^2 - q^2) / 2 / k / p
    if (abs(abs(legendre_x) - 1) < 1e-12)
        legendre_x = sign(legendre_x) * 1
    end
    # convention changed, now interaction_dynamic already included the V_Bare
    return q * Pl(legendre_x, channel) * interaction_dynamic(q, n, param, int_type, spin_state)
end

@inline function kernel0_integrand(k, p, q, channel, param, spin_state)
    legendre_x = (k^2 + p^2 - q^2) / 2 / k / p
    if (abs(abs(legendre_x) - 1) < 1e-12)
        legendre_x = sign(legendre_x) * 1
    end
    @assert -1 <= legendre_x <= 1 "k=$k,p=$p,q=$q"

    return q * Pl(legendre_x, channel) * interaction_instant(q, param, spin_state)
end

@inline function kernel_integrand2d(k, p, θ, n, channel, param, int_type, spin_state)
    x = cos(θ)
    q2 = k^2 + p^2 - 2k * p * x
    if x == 1 || q2 <= 0
        q2 = (k - p)^2 + k * p * θ^2
    end
    return cos(channel * θ) * interaction_dynamic(√q2, n, param, int_type, spin_state)
end

@inline function kernel0_integrand2d(k, p, θ, channel, param, spin_state)
    x = cos(θ)
    q2 = k^2 + p^2 - 2k * p * x
    if x == 1 || q2 <= 0
        q2 = (k - p)^2 + k * p * θ^2
        # println("$k, $p, $q2, $((k - p)^2), $x, $θ")
    end
    return cos(channel * θ) * interaction_instant(√q2, param, spin_state)
end

function helper_function(y::Float64, n::Int, W, param; Nk::Int = 40, minK::Float64 = 1e-12 * param.kF, order::Int = 6)
    # return the helper function
    @unpack kF, β = param
    # generate a new grid for every calculation
    kgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, y], [0.0, min(y, 2kF)], Nk, minK, order)

    integrand = zeros(Float64, kgrid.size)
    for (ki, k) in enumerate(kgrid)
        integrand[ki] = k^n * W(k)
    end

    H = Interp.integrate1D(integrand, kgrid)
    return H
end

function helper_function_grid(ygrid, intgrid, n::Int, W, param)
    # return the helper function
    @unpack kF, β, e0, ϵ0 = param
    # generate a new grid for every calculation
    #kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, grid[end]], [0.0,min(grid[end],2kF)], Nk, minK, order)
    kgrid = intgrid
    grid = ygrid
    helper = zeros(Float64, length(grid))

    integrand = zeros(Float64, kgrid.size)
    for (ki, k) in enumerate(kgrid)
        integrand[ki] = k^(n-2) * W(k) * e0^2 / ϵ0
    end

    for i in 1:length(grid)
        if i == 1
            x1, x2, hprev = 0.0, grid[1], 0.0
        else
            x1, x2, hprev = grid[i-1], grid[i], helper[i-1]
        end
        helper[i] = Interp.integrate1D(integrand, kgrid, [x1, x2]) + hprev
        # helper[i] = Interp.integrate1D(integrand, kgrid, [EPS,grid[i]])
        # @assert isfinite(helper[i]) "fail at $(grid[i])"
    end

    return helper
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

    function DCKernel(int_type, spin_state, channel, param, kgrid, qgrids, dlrGrid, kernel_bare, kernel)
        return new(int_type, spin_state, channel, param, kgrid, qgrids, dlrGrid, kernel_bare, kernel)
    end

end

function DCKernel_2d(param, Euv, rtol, Nk, maxK, minK, order, int_type, channel, spin_state = :auto)
    @unpack kF, β = param

    if spin_state == :sigma
        # for self-energy, always use ℓ=0
        channel = 0
    elseif spin_state == :auto
        # automatically assign spin_state, triplet for even, singlet for odd channel
        spin_state = (channel % 2 == 0) ? (:triplet) : (:singlet)
    end

    bdlr = DLRGrid(Euv, β, rtol, false, :ph)
    kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order)
    #println(kgrid.grid)
    qgrids = [CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [k, kF], Nk, minK, order) for k in kgrid.grid]
    qgridmax = maximum([qg.size for qg in qgrids])
    # θgrid = SimpleGrid.GaussLegendre{Float64}([0, 2π], 100)
    # θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], Nk, 1e-7, order)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], Nk, 1e-7 * kF, order)
    # println(θgrid.size)
    # println(θgrid.grid)

    kernel_bare = zeros(Float64, (length(kgrid.grid), (qgridmax)))
    kernel = zeros(Float64, (length(kgrid.grid), (qgridmax), length(bdlr.n)))
    data = zeros(Float64, length(θgrid.grid))
    for (ki, k) in enumerate(kgrid.grid)
        for (pi, p) in enumerate(qgrids[ki].grid)
            for (θi, θ) in enumerate(θgrid.grid)
                data[θi] = kernel0_integrand2d(k, p, θ, channel, param, spin_state)
            end
            kernel_bare[ki, pi] = Interp.integrate1D(data, θgrid) * 2
            # println("$ki,$pi,  $(kernel_bare[ki,pi])")
            @assert isfinite(kernel_bare[ki, pi]) "fail kernel_bare at $ki,$pi, ($k, $p) with $(kernel_bare[ki,pi]) \n $data"
            for (ni, n) in enumerate(bdlr.n)
                for (θi, θ) in enumerate(θgrid.grid)
                    data[θi] = kernel_integrand2d(k, p, θ, n, channel, param, int_type, spin_state)
                end
                # data = [kernel_integrand2d(k, p, θ, n, channel, param, int_type, spin_state) for θ in θgrid.grid]
                kernel[ki, pi, ni] = Interp.integrate1D(data, θgrid) * 2
                # ni == 1 && println("$ki,$pi, $n  $(kernel[ki, pi, ni])")
                @assert isfinite(kernel[ki, pi, ni]) "fail kernel at $ki,$pi,$ni, with $(kernel[ki,pi,ni])"
            end
        end
    end

    return DCKernel(int_type, spin_state, channel, param, kgrid, qgrids, bdlr, kernel_bare, kernel)
end

function DCKernel_2d(param; Euv = param.EF * 100, rtol = 1e-10, Nk = 8, maxK = param.kF * 10, minK = param.kF * 1e-7, order = 4, int_type = :rpa, channel = 0, spin_state = :auto)
    return DCKernel_2d(param, Euv, rtol, Nk, maxK, minK, order, int_type, channel, spin_state)
end

function DCKernel_old(param, Euv, rtol, Nk, maxK, minK, order, int_type, channel, spin_state = :auto)
    @unpack kF, β = param

    if spin_state == :sigma
        # for self-energy, always use ℓ=0
        channel = 0
    elseif spin_state == :auto
        # automatically assign spin_state, triplet for even, singlet for odd channel
        spin_state = (channel % 2 == 0) ? (:triplet) : (:singlet)
    end

    bdlr = DLRGrid(Euv, β, rtol, false, :ph)
    kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order)
    #println(kgrid.grid)
    qgrids = [CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [k, kF], Nk, minK, order) for k in kgrid.grid]
    qgridmax = maximum([qg.size for qg in qgrids])
    #println(qgridmax)

    kernel_bare = zeros(Float64, (length(kgrid.grid), (qgridmax)))
    kernel = zeros(Float64, (length(kgrid.grid), (qgridmax), length(bdlr.n)))

    int_grid_base = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.1 * maxK], [0.0, 2kF], 2Nk, 0.01minK, 2)
    for (ki, k) in enumerate(kgrid.grid)
        for (pi, p) in enumerate(qgrids[ki].grid)
            if abs(k - p) > EPS

                kmp = abs(k - p) < EPS ? EPS : abs(k - p)
                kpp = k + p
                im, ip = floor(int_grid_base, kmp), floor(int_grid_base, kpp)
                int_panel = Float64[]

                push!(int_panel, kmp)
                if im < ip
                    for i in im+1:ip
                        push!(int_panel, int_grid_base[i])
                    end
                end
                push!(int_panel, kpp)

                int_panel = SimpleGrid.Arbitrary{Float64}(int_panel)
                SubGridType = SimpleGrid.GaussLegendre{Float64}
                subgrids = subgrids = Vector{SubGridType}([])
                for i in 1:int_panel.size-1
                    _bound = [int_panel[i], int_panel[i+1]]
                    push!(subgrids, SubGridType(_bound, order))
                end

                int_grid = CompositeGrid.Composite{Float64,SimpleGrid.Arbitrary{Float64},SubGridType}(int_panel, subgrids)

                data = [kernel0_integrand(k, p, q, channel, param, spin_state) for q in int_grid.grid]
                kernel_bare[ki, pi] = Interp.integrate1D(data, int_grid)

                for (ni, n) in enumerate(bdlr.n)
                    data = [kernel_integrand(k, p, q, n, channel, param, int_type, spin_state) for q in int_grid.grid]
                    kernel[ki, pi, ni] = Interp.integrate1D(data, int_grid)
                    @assert isfinite(kernel[ki, pi, ni]) "fail kernel at $ki,$pi,$ni, with $(kernel[ki,pi,ni])"
                end

            else
                kernel_bare[ki, pi] = 0
                for (ni, n) in enumerate(bdlr.n)
                    kernel[ki, pi, ni] = 0
                end
            end
        end
    end

    return DCKernel(int_type, spin_state, channel, param, kgrid, qgrids, bdlr, kernel_bare, kernel)
end

function DCKernel_old(param; Euv = param.EF * 100, rtol = 1e-10, Nk = 8, maxK = param.kF * 10, minK = param.kF * 1e-7, order = 4, int_type = :rpa, channel = 0, spin_state = :auto)
    return DCKernel_old(param, Euv, rtol, Nk, maxK, minK, order, int_type, channel, spin_state)
end

function DCKernel0(param; Euv = param.EF * 100, rtol = 1e-10, Nk = 8, maxK = param.kF * 10, minK = param.kF * 1e-7, order = 4, int_type = :rpa, spin_state = :auto)
    return DCKernel0(param, Euv, rtol, Nk, maxK, minK, order, int_type, spin_state)
end

function DCKernel0(param, Euv, rtol, Nk, maxK, minK, order, int_type, spin_state = :auto)
    # use helper function
    @unpack kF, β = param
    channel = 0

    if spin_state == :sigma
        # for self-energy, always use ℓ=0
        channel = 0
    elseif spin_state == :auto
        # automatically assign spin_state, triplet for even, singlet for odd channel
        spin_state = (channel % 2 == 0) ? (:triplet) : (:singlet)
    end

    bdlr = DLRGrid(Euv, β, rtol, false, :ph)
    kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order)
    #println(kgrid.grid)
    qgrids = [CompositeGrid.LogDensedGrid(:gauss, [0.0, maxK], [k, kF], Nk, minK, order) for k in kgrid.grid]
    qgridmax = maximum([qg.size for qg in qgrids])
    #println(qgridmax)

    kernel_bare = zeros(Float64, (length(kgrid.grid), (qgridmax)))
    kernel = zeros(Float64, (length(kgrid.grid), (qgridmax), length(bdlr.n)))

    helper_grid = CompositeGrid.LogDensedGrid(:cheb, [0.0, 2.1 * maxK], [0.0, 2kF], 2Nk, 0.01minK, 2order)
    intgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, helper_grid[end]], [0.0, 2kF], 2Nk, 0.01minK, 2order)

    # dynamic
    for (ni, n) in enumerate(bdlr.n)
        # helper = zeros(Float64, helper_grid.size)
        # for (yi, y) in enumerate(helper_grid)
        #     helper[yi] = helper_function(y, 1, u->interaction_dynamic(u,n,param,int_type,spin_state),param)
        # end
        helper = helper_function_grid(helper_grid, intgrid, 1, u -> interaction_dynamic(u, n, param, int_type, spin_state), param)
        for (ki, k) in enumerate(kgrid.grid)
            for (pi, p) in enumerate(qgrids[ki].grid)
                Hp, Hm = Interp.interp1D(helper, helper_grid, k + p), Interp.interp1D(helper, helper_grid, abs(k - p))
                kernel[ki, pi, ni] = (Hp - Hm)
            end
        end
    end

    # instant
    # helper = zeros(Float64, helper_grid.size)
    # for (yi, y) in enumerate(helper_grid)
    #     helper[yi] = helper_function(y, 1, u->interaction_instant(u,param,spin_state),param)
    # end

    # helper = helper_function_grid(helper_grid, intgrid, 1, u -> interaction_instant(u, param, spin_state), param)
    helper = helper_function_grid(helper_grid, intgrid, 1, u -> 1.0, param)
    for (ki, k) in enumerate(kgrid.grid)
        for (pi, p) in enumerate(qgrids[ki].grid)
            Hp, Hm = Interp.interp1D(helper, helper_grid, k + p), Interp.interp1D(helper, helper_grid, abs(k - p))
            kernel_bare[ki, pi] = (Hp - Hm)
        end
    end

    return DCKernel(int_type, spin_state, channel, param, kgrid, qgrids, bdlr, kernel_bare, kernel)
end

end
