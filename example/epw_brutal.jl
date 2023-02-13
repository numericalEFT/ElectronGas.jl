
# use brutal force way to compute Tc from EPW data

include("epw_io.jl")

fermi_mat(n, β, N) = π * (2(n - (N ÷ 2) - 1) + 1) / β

function s_matrix(N, wsph, a2f_iso, β, muc=0.07)
    smat = zeros(Float64, (N, N))
    for i in 1:N
        for j in 1:N
            wi, wj = fermi_mat(i, β, N), fermi_mat(j, β, N)
            λ = lambdar_iso(wi - wj, wsph, a2f_iso)
            smat[i, j] = λ - muc
            if i == j
                for k in 1:N
                    wk = fermi_mat(k, β, N)
                    smat[i, j] = smat[i, j] - lambdar_iso(wi - wk, wsph, a2f_iso) * sign(wi * wk)
                end
            end
            smat[i, j] = smat[i, j] / abs(wi) * π / β
        end
    end
    return smat
end

sqnorm(x) = sqrt(sum(x .^ 2))

function power_method(smat; conv=1e-6, Nmax=200, shift=1.0)
    N = size(smat)[1]
    x, y = zeros(Float64, N), ones(Float64, N)

    a = 0.0
    diff = 1.0
    n = 1
    while (diff > conv && n < Nmax)
        norm = sqnorm(y)
        x = y ./ norm
        # y = y .* 0
        # for i in 1:N
        #     for j in 1:N
        #         y[i] = y[i] + smat[i, j] * x[j]
        #     end
        # end
        y = smat * x
        y = y .+ shift .* x
        diff = abs(a - norm)
        a = norm
        # a = 0.0
        # for i in 1:N
        #     a = a + x[i] * y[i]
        # end
        # x = y .- a .* x
        # norm = sqnorm(x)
        # diff = norm
        # println(n, ",", a - shift, ",", y[N÷2], ",", y[1])

        n = n + 1
    end
    return a - shift
end

function compute_λ(T, Ec, wsph, a2f_iso)
    β = 1 / T
    N = floor(Int, 2Ec / T)
    smat = s_matrix(N, wsph, a2f_iso, β)
    λ = power_method(smat)
    println("λ=", λ, ", at T=", 1.160451812e4 / β)
    return λ
end

@testset "QE_BRUTAL" begin
    prefix = "pb"
    # dir = "~/File/Research/Quantum-Espresso/EPW/Thu.6.Margine/exercise1/epw/"
    dir = "./run/epw/"

    wsph, a2f_iso = read_a2f(prefix; dir=dir)

    compute_λ(0.00044, 0.1, wsph, a2f_iso)
    compute_λ(0.00045, 0.1, wsph, a2f_iso)
    compute_λ(0.00046, 0.1, wsph, a2f_iso)
end