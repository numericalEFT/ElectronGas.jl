using ElectronGas, JLD2
using ElectronGas.GreenFunc
using ElectronGas.CompositeGrids

function load_bs(; uid0=3020000, n=84)
    Bs = []
    uids = []
    param = nothing
    for duid in 1:n
        uid = uid0 + duid
        fname = "run/data/calc_b/PCFresumdlr_$uid.jld2"
        try
            f = jldopen(fname, "r")
            B = f["B"]
            param = f["param"]
            close(f)
            push!(Bs, B)
            push!(uids, uid)
        catch e
            println(e)
            println("uid=$uid")
            uid = uids[end]
            fname = "run/data/calc_b/PCFresumdlr_$uid.jld2"
            f = jldopen(fname, "r")
            B = f["B"]
            param = f["param"]
            close(f)
            push!(Bs, B)
            push!(uids, uid)
        end
    end
    B = similar(Bs[1])
    for i in 1:length(Bs)
        B.data[i, :] .= Bs[i].data[i, :]
    end
    return param, B
end

function load_A(fname)
    f = jldopen(fname, "r")
    A = f["A"]
    param = f["param"]
    return param, A
end

function interp_A(param, A, B; Ec=100param.EF)
    # interp A to wgrid of B
    Adlr = to_dlr(A)
    nec = floor(Int, Ec / 2 / π * param.β + 0.5)
    Adense = dlr_to_imfreq(Adlr, [i for i in 0:nec])
    densegrid = SimpleG.Arbitrary([Adense.mesh[1][i] for i in 1:length(Adense)])
    wgrid = B.mesh[1]
    newA = MeshArray(wgrid; dtype=Float64)
    for i in 1:length(newA)
        w = newA.mesh[1][i]
        newA.data[i] = Interp.interp1D(Adense.data, densegrid, w)
    end
    return newA
end

function save_AB(fname, param, A, B)
    f = jldopen(fname, "w")
    f["A"] = A
    f["param"] = param
    f["B"] = B
    close(f)
end

# param, A = load_A("run/data/PCFresumdlr_3000022.jld2")
# param, A = load_A("run/data/Adlr_ko3_beta6400_lam6.jld2")
param, A = load_A("run/data/calc_b/PCFresumdlr_30600000.jld2")
param, B = load_bs(uid0=30600000)

# println(B.data)

# newA = interp_A(param, A, B)
newA = A
save_AB("run/data/PCFresumrs3_3066022.jld2", param, newA, B)
