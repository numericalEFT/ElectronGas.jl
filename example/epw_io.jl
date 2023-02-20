# this script test loading data from epw

module QE_EPW

using DelimitedFiles

using ElectronGas
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids
using ElectronGas.Parameters
using ElectronGas.Lehmann

using ElectronGas.Parameter
using ElectronGas.BSeq

export read_a2f, lambdar_iso, lambda_wrapped, G2ωwrapped, linearResponse, measure_R

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
                wsph[line] = wsph[line] / 1000.0 # EPW did this, meV to eV
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
        green_dyn[ni] = -π * lambdar_iso(wn_mesh[ni], wsph, a2f_iso)
    end
    green_ins[1] = muc * π
    # green_ins[1] = g
    return green_dyn, green_ins
end

function G2ωwrapped(Euv, rtol, param)
    @unpack β, Ωcut = param
    wn_mesh = GreenFunc.ImFreq(β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    # print("max frequency $(wn_mesh[end]*β)\n")
    green = GreenFunc.MeshArray(wn_mesh; dtype=Float64)
    for ind in eachindex(green)
        ni = ind[1]
        n = wn_mesh.grid[ni]
        wn = wn_mesh[ni]
        green[ind] = 1 / abs(wn)
        # green[ind] = π / abs(wn) # EPW convention
    end
    return green
end

function Δω_correction(G2, wsph, a2f_iso, nmax)
    correction = G2 * 0.0
    for ind in eachindex(correction)
        ni = ind[1]
        n = G2.mesh[1].grid[ni]
        wn = G2.mesh[1][ni]
        for np in (1-nmax):nmax
            wnp = (2np + 1) * π / G2.mesh[1].β
            correction[ni] -= sign(wnp) * sign(wn) * lambdar_iso(wn - wnp, wsph, a2f_iso) / abs(2n + 1)
        end
        if n == nmax
            correction[ni] = 0.0
        end
    end
    return correction
end

function linearResponse(param, wsph, a2f_iso; Ecut=100000 * param.EF,
    Euv=Ecut / param.β, rtol=1e-10, atol=1e-10, α=0.8, verbose=false, Ntherm=30, muc=0.1)
    @unpack β = param
    @time kernel_freq, kernel_ins = lambda_wrapped(Euv, rtol, param, wsph, a2f_iso; muc=muc)
    G2 = G2ωwrapped(Euv, rtol, param)
    fdlr = Lehmann.DLRGrid(Euv, β, rtol, FERMION, :pha)
    # bdlr = kernel_freq.mesh[1].dlrGrid
    @assert G2.mesh[1].grid == fdlr.n
    # prepare kernel, interpolate into τ-space with fdlr.τ
    # kernel =  kernel_freq |> to_dlr |>to_imtime
    kernel = dlr_to_imtime(to_dlr(kernel_freq), fdlr.τ)

    # delta_correction corresponds to additional term in s matrix
    # delta_correction = G2 * 0.0
    nmax = G2.mesh[1].grid[end]
    delta_correction = Δω_correction(G2, wsph, a2f_iso, nmax)
    # println(delta_correction.data)

    # calculate F, R by Bethe-Slapter iteration.
    lamu, F_freq, R_imt, R_ins = BSeq.BSeq_solver(param, G2, kernel, kernel_ins, Euv;
        rtol=rtol, α=α, atol=atol, verbose=verbose, Ntherm=Ntherm, delta_correction=delta_correction)
    # println("1/R₀ = $lamu, at T = $(157887.663431925/β)")
    println("1/R₀ = $lamu, at T = $(1.160451812e4/β)")

    R_freq = R_imt |> to_dlr |> to_imfreq
    # println(view(R_freq, :, kF_label))
    return lamu, R_freq, F_freq
end

function measure_R(betas, prefix;
    dir="./run/epw/", Ecut=160000, rtol=1e10,
    kwargs...)

    wsph, a2f_iso = read_a2f(prefix; dir=dir)
    lamus = betas .* 0.0

    for i in 1:length(betas)
        param = QE_EPW.Parameter.defaultUnit(1 / betas[i], 1.0)
        lambda_dyn, lambda_ins = lambda_wrapped(Ecut / betas[i], rtol, param, wsph, a2f_iso)
        lamu, R_freq, F_freq = linearResponse(param, wsph, a2f_iso; Ecut=Ecut, kwargs...)
        lamus[i] = lamu

        # measure chi
        F_dlr = F_freq |> to_dlr
        F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
        chi = F_ins[1]
        # write 
        data = [betas[i] 1 / chi lamu 0 1.0]

        dir = "./run/"
        fname = "epw_chi_rs$(1.0)_l$(0).txt"
        # fname = "gap$(dim)D_chi_rs$(rs)_l$(channel).txt"
        open(dir * fname, "a+") do io
            writedlm(io, data, ' ')
        end
    end

    return lamus
end

end

using Test
using Lehmann
using LinearAlgebra
using .QE_EPW

const ev2Kelvin = 1.160451812e4
function linreg(X, Y)
    hcat(fill!(similar(X), 1), X) \ Y
end

if abspath(PROGRAM_FILE) == @__FILE__
    @testset "QE_EPW" begin

        # # test read_a2f

        # prefix = "pb"
        # # dir = "~/File/Research/Quantum-Espresso/EPW/Thu.6.Margine/exercise1/epw/"
        # dir = "./run/epw/"

        # wsph, a2f_iso = read_a2f(prefix; dir=dir)
        # # println(wsph)
        # # println(a2f_iso)

        # # test lambdar_iso

        # println("lambdar_iso(1.0)=", lambdar_iso(1.0, wsph, a2f_iso))

        # # test DLR
        # Euv = 100
        # β = 1000
        # rtol = 1e-10
        # bdlr = DLRGrid(Euv, β, rtol, false, :ph)

        # λ = [lambdar_iso(bdlr.ωn[i], wsph, a2f_iso) for i in 1:length(bdlr)]

        # # println(bdlr.ωn)
        # # println(λ)

        # # test lambda_wrapped
        # param = QE_EPW.Parameter.defaultUnit(0.0007, 1.0)
        # lambda_dyn, lambda_ins = lambda_wrapped(Euv, rtol, param, wsph, a2f_iso)

        # println(lambda_dyn.mesh[1])
        # println(lambda_dyn.data)

        # # test linear solver
        # lamu, R_freq, F_freq = linearResponse(param, wsph, a2f_iso)

        # param = QE_EPW.Parameter.defaultUnit(0.0006, 1.0)
        # lambda_dyn, lambda_ins = lambda_wrapped(Euv, rtol, param, wsph, a2f_iso)
        # lamu, R_freq, F_freq = linearResponse(param, wsph, a2f_iso)

        num = 9
        # betas = [6.25 * sqrt(2)^i for i in LinRange(0, num - 1, num)]
        # betas = [1 / 0.0005, 1 / 0.00055, 1 / 0.0006]
        betas = [1 / (5.5 * 2^((i - 1) / num)) for i in 1:num]
        betas = betas .* ev2Kelvin
        lamus = measure_R(betas, "pb"; muc=0.1)
        println(betas)
        println(lamus)
        lnbetas = log10.(betas)
        b1, k1 = linreg(lnbetas, lamus)
        println("k=$k1, b=$b1, Tc=$(10^(b1/k1)*ev2Kelvin)")

    end
end