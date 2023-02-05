# this script test loading data from epw

module QE_EPW

export read_a2f, lambdar_iso

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

    println(bdlr.ωn)
    println(λ)
end