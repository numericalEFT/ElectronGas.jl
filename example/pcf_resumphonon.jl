using ElectronGas
using JLD2
using ElectronGas.CompositeGrids
using ElectronGas.Parameters
using ElectronGas.GreenFunc

function load_B(fname)
    f = jldopen(fname, "r")
    param = f["param"]
    B = f["B"]
    return param, B
end

function save_B(fname, param, B)
    f = jldopen(fname, "w")
    f["param"] = param
    f["B"] = B
    close(f)
end

function avg_phonon(n, param)
    # compute phonon averaged on fermi surface
    xgrid = SimpleG.GaussLegendre{Float64}([-1.0, 1.0], 20)
    data = zeros(Float64, length(xgrid))
    for i in 1:length(xgrid)
        data[i] = Interaction.phonon(n, param.kF * sqrt(2 - 2 * xgrid[i]), param)[1]
    end
    return Interp.integrate1D(data, xgrid)
end

function add_phonon(B, param)
    # add phonon to B
    wgrid = B.mesh[1]
    for (wi, w) in enumerate(wgrid)
        for (vi, v) in enumerate(wgrid)
            dw = abs(w - v)
            n = dw / 2 / π * param.β
            B.data[wi, vi] -= avg_phonon(n, param)
            dw = abs(w + v)
            n = dw / 2 / π * param.β
            B.data[wi, vi] -= avg_phonon(n, param)
        end
    end
    return B
end

fname = "run/data/PCFresumrs3_3055022.jld2"
param, B = load_B(fname)
param = Parameter.Para(param; eph=0.4 * 4 * π^2, ω_D=0.005 * param.EF)
Bph = add_phonon(B, param)

savename = "run/data/Bsmooth_ko3ph4_beta6400_lam6.jld2"
save_B(savename, param, B)