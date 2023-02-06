# this script test loading data from epw

module QE_EPW

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids
using ElectronGas.Parameters
using ElectronGas.Lehmann

using ElectronGas.Parameter
using ElectronGas.BSeq

export read_a2f, lambdar_iso, lambda_wrapped, G2ωwrapped, linearResponse

function read_a2f(prefix; dir=nothing, suffix=nothing, nqstep=500)
    # read_a2f in EPW/src/io_eliashberg.f90
    if isnothing(dir)
        dir = "./"
    end
    if isnothing(suffix)
        suffix = ".a2f"
    end
    fname = dir * prefix * suffix
    println("loading $fname")

    wsph = zeros(Float64, nqstep)
    a2f_iso = zeros(Float64, nqstep)

    open(fname) do f
        line = 0
        while !eof(f)
            s = readline(f)
            if line in 1:nqstep
                ssplit = split(s)
                wsph[line] = parse(Float64, ssplit[1])
                a2f_iso[line] = parse(Float64, ssplit[2])
                wsph[line] = wsph[line] / 1000.0 # EPW did this
            end
            line += 1
        end
    end

    return wsph, a2f_iso
end

function lambdar_iso(omega, wsph, a2f_iso)
    lambda_eph = 0.0
    nqstep = length(wsph)
    dwsph = wsph[end] / nqstep

    for iwph in 1:nqstep
        lambda_eph = lambda_eph + wsph[iwph] * a2f_iso[iwph] / (wsph[iwph] * wsph[iwph] + omega^2)
    end

    lambda_eph = 2.0 * lambda_eph * dwsph

    return lambda_eph
end

function lambda_wrapped(Euv, rtol, param, wsph, a2f_iso; muc=0.1)
    @unpack β, g = param
    wn_mesh = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn = GreenFunc.MeshArray(wn_mesh; dtype=Float64)
    green_ins = GreenFunc.MeshArray([0.0,]; dtype=Float64)
    for (ni, n) in enumerate(wn_mesh.grid)
        # green_dyn[ni]=gamma(n,param) 
        green_dyn[ni] = -lambdar_iso(wn_mesh[ni], wsph, a2f_iso)
    end
    green_ins[1] = muc
    # green_ins[1] = g
    return green_dyn, green_ins
end

function G2ωwrapped(Euv, rtol, param)
    @unpack β, Ωcut = param
    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    print("max frequency $(wn_mesh[end]*β)\n")
    green = GreenFunc.MeshArray(wn_mesh; dtype=Float64)
    for ind in eachindex(green)
        ni = ind[1]
        n = wn_mesh.grid[ni]
        wn = wn_mesh[ni]
        green[ind] = 1.0 / abs(wn)
    end
    return green
end

function linearResponse(param, wsph, a2f_iso; Euv=100 * param.EF, rtol=1e-10, atol=1e-10, α=0.7, verbose=false, Ntherm=30, muc=0.1)
    @unpack β = param
    @time kernel_freq, kernel_ins = lambda_wrapped(Euv, rtol, param, wsph, a2f_iso; muc=muc)
    G2 = G2ωwrapped(Euv, rtol, param)
    fdlr = Lehmann.DLRGrid(Euv, β, rtol, FERMION, :pha)
    # bdlr = kernel_freq.mesh[1].dlrGrid
    @assert G2.mesh[1].grid == fdlr.n
    # prepare kernel, interpolate into τ-space with fdlr.τ
    # kernel =  kernel_freq |> to_dlr |>to_imtime
    kernel = dlr_to_imtime(to_dlr(kernel_freq), fdlr.τ)

    # calculate F, R by Bethe-Slapter iteration.
    lamu, F_freq, R_imt, R_ins = BSeq.BSeq_solver(param, G2, kernel, kernel_ins, Euv;
        rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm)
    println("1/R₀ = $lamu, at T = $(157887.663431925/β)")

    R_freq = R_imt |> to_dlr |> to_imfreq
    # println(view(R_freq, :, kF_label))
    return lamu, R_freq, F_freq
end

end

using Test
using Lehmann
using .QE_EPW

@testset "QE_EPW" begin

    # test read_a2f

    prefix = "pb"
    # dir = "~/File/Research/Quantum-Espresso/EPW/Thu.6.Margine/exercise1/epw/"
    dir = "./run/epw/"

    wsph, a2f_iso = read_a2f(prefix; dir=dir)
    # println(wsph)
    # println(a2f_iso)

    # test lambdar_iso

    println("lambdar_iso(1.0)=", lambdar_iso(1.0, wsph, a2f_iso))

    # test DLR
    Euv = 100
    β = 1000
    rtol = 1e-10
    bdlr = DLRGrid(Euv, β, rtol, false, :ph)

    λ = [lambdar_iso(bdlr.ωn[i], wsph, a2f_iso) for i in 1:length(bdlr)]

    # println(bdlr.ωn)
    # println(λ)

    # test lambda_wrapped
    param = QE_EPW.Parameter.defaultUnit(0.00021, 1.0)
    lambda_dyn, lambda_ins = lambda_wrapped(Euv, rtol, param, wsph, a2f_iso)

    println(lambda_dyn.mesh[1])
    println(lambda_dyn.data)

    # test linear solver
    lamu, R_freq, F_freq = linearResponse(param, wsph, a2f_iso)
end