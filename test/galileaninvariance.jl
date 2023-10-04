# using ElectronGas
using Parameters
using GreenFunc
using CompositeGrids
# print(dirname(pwd()))
include("./ElectronGas.jl/src/convention.jl")
include("./ElectronGas.jl/src/parameter.jl")
include("./ElectronGas.jl/src/polarization.jl")
include("./ElectronGas.jl/src/interaction.jl")
include("./ElectronGas.jl/src/twopoint.jl")
using Plots

para = Parameter.rydbergUnit(0.01, 1, 3)
Ωn_mesh_ph = GreenFunc.ImFreq(para.β, BOSON; Euv=100.0, symmetry=:ph)
sgrid = 0.0

δladderarray_n = GreenFunc.MeshArray(Ωn_mesh_ph, sgrid; dtype=ComplexF64)
δTarray_n2 = GreenFunc.MeshArray(Ωn_mesh_ph, sgrid; dtype=ComplexF64)

for (ki, k) in enumerate(sgrid)
    for (ni, n) in enumerate(Ωn_mesh_ph.grid)
        δTarray_n2[ni, ki] = Interaction.Tmatrix(k, n, para) - Interaction.T0matrix(k, n, para)
    end
end




# Many Body Particle-Particle
mbladderarray = Polarization.Ladder0_FiniteTemp(0.0, Ωn_mesh_ph.grid, para)

# Vacuum Particle-Particle
vacladderarray = Polarization.Ladder0_FreeElectron(0.0,  Ωn_mesh_ph.grid, para)

for (ki, k) in enumerate(sgrid)
    for (ni, n) in enumerate(Ωn_mesh_ph.grid)
        δladderarray_n[ni, ki] = mbladderarray[ni,ki] - vacladderarray[ni,ki] 
    end
end

# for (ki, k) in enumerate(sgrid)
#     for (ni, n) in enumerate(Ωn_mesh_ph.grid)
#         δTarray_n[ni, ki] = Interaction.Tmatrix(k, n, para) - Interaction.T0matrix(k, n, para)
#     end
# end

δladderarray_τ = δladderarray_n |> to_dlr |> to_imtime
# δTarray_τ = δTarray_n |> to_dlr |> to_imtime

# Analytic expression for T0(k, τ) #

function T0_τ(q::Float64, τ::Float64, param)

    me, as, μ, β  = param.me, param.as, param.μ, param.β
    B = me^(3/2) / (4 * π)
    # Generating a log densed composite grid with LogDensedGrid()
    ygrid =  CompositeGrid.LogDensedGrid(:gauss,[0.0,50.],[0.0],5,0.005,5)

    # Functions to be integrated 
    f(y) = ( exp(- y^2 ) * y^2 ) / (τ / (me * as) + y^2)
    g(y) = ( exp(- y^2 ) * y^2 ) / (τ / (me * as) + y^2) * 1. / ( exp( β * ( y^2/ τ + (q^2 / (4me) - 2μ)) ) - 1)

    # Generate data for integration
    fdata = [f(y) for (yi, y) in enumerate(ygrid.grid)]
    gdata = [g(y) for (yi, y) in enumerate(ygrid.grid)]

    fint_result = Interp.integrate1D(fdata, ygrid)
    gint_result = Interp.integrate1D(gdata, ygrid)

    return exp(-(q^2 / (4 * me) - 2μ) * τ) * (2/( B* π * sqrt(τ) ) * ( fint_result + gint_result ) - exp( τ / (me*as^2) ) / ( 1 - exp( β* ( 1 / (me*as^2) - (q^2 / (4me) - 2μ)) ) ))
    # return  - exp( τ / (me*as^2) ) / ( 1 - exp( β* ( 1 / (me*as^2) - (q^2 / (4me) - 2μ)) ) )

end

function Σ_τ(q::Float64, param)
    # Output self energy for the mesh imaginary values

    kgrid = Polarization.finitetemp_kgrid_ladder(param.μ, param.me, q, param.kF, 20, 1e-6, 10)
    δTarray_n = GreenFunc.MeshArray(Ωn_mesh_ph, kgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(Ωn_mesh_ph.grid)
            δTarray_n[ni, ki] = Interaction.Tmatrix(k+q, n, param) - Interaction.T0matrix(k+q, n, param)
        end
    end

    δTarray_τ = δTarray_n |> to_dlr |> to_imtime
    integrand = GreenFunc.MeshArray(δTarray_τ, kgrid; dtype=ComplexF64)
    result = GreenFunc.MeshArray(δTarray_τ.mesh[1]; dtype=ComplexF64)

    for (ti, t) in enumerate(δTarray_τ.mesh[1][:])
        for (ki, k) in enumerate(kgrid.grid)
            integrand[ti, ki] = (T0_τ(k + q, t, param) + δTarray_τ[ti]) * exp(t * ( k^2 / (2*param.me)  - param.μ)) * 1 / (exp( param.β * (k^2 / (2* param.me) - param.μ)) + 1)
        end
        result[ti] = Interp.integrate1D(integrand[ti, :], kgrid)
    end
    return result
end

qarray = LinRange(0,2*para.kF, 5)
#test_τ = #[Σ_τ(qi, para) for qi in qarray]
#test_ω = test_τ |> to_dlr |> to_imfreq
#plot(test_ω.mesh[1][:], imag(test_ω))

#G0_τ = GreenFunc.MeshArray(test_τ.mesh[1]; dtype=ComplexF64)

#for (ti, t) in enumerate(G0_τ.mesh[1][:])
#    G0_τ[ti] = T0_τ(0.0, t, para) #exp(t * ( 0.0^2 / (2*para.me)  - para.μ)) * 1 / (exp( para.β * (0.0^2 / (2* para.me) - para.μ)) + 1)

#end
#G0_ω = G0_τ |> to_dlr |> to_imfreq