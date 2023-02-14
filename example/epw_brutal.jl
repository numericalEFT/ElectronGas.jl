
# use brutal force way to compute Tc from EPW data

include("epw_io.jl")

const ev2Kelvin = 1.160451812e4

fermi_mat(n, β, N) = π * (2(n - (N ÷ 2) - 1) + 1) / β

function s_matrix(N, wsph, a2f_iso, β, muc=0.1)
    # println("$N, $(fermi_mat(1, β, N)), $(fermi_mat(N÷2, β, N)), $(fermi_mat(N, β, N)) ")
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

function precursory_cooper_flow(smat;
    source=ones(Float64, size(smat)[1]), α=0.8,
    conv=1e-6, Nmax=1000)

    N = size(smat)[1]

    invR0 = 0.0
    diff = 1.0
    n = 1
    xsum = source
    x = xsum .* (1 - α)

    converge = false

    while (!converge && n < Nmax)
        y = smat * x
        xsum = xsum .* α .+ (source + y)
        x = xsum .* (1 - α)
        # diff = abs(1 / x[N÷2] - invR0)
        converge = isapprox(1 / x[N÷2], invR0, rtol=conv, atol=1e-10)
        invR0 = 1 / x[N÷2]
        n = n + 1
    end

    return invR0
end

function compute_λ(T, Ec, wsph, a2f_iso)
    β = 1 / T
    N = floor(Int, Ec / π / T / 2) * 2
    smat = s_matrix(N, wsph, a2f_iso, β)
    λ = power_method(smat)
    println("λ=", λ, ", at T=", 1.160451812e4 / β)
    println("λ from LA: $(eigmax(smat))")
    return λ
end

function compute_invR0(T, Ec, wsph, a2f_iso)
    β = 1 / T
    N = floor(Int, Ec / π / T / 2) * 2
    smat = s_matrix(N, wsph, a2f_iso, β)
    invR0 = precursory_cooper_flow(smat)
    println("1/R0=", invR0, ", at T=", ev2Kelvin / β)
    return invR0
end

function linreg(X, Y)
    hcat(fill!(similar(X), 1), X) \ Y
end

@testset "QE_BRUTAL" begin
    prefix = "pb"
    # dir = "~/File/Research/Quantum-Espresso/EPW/Thu.6.Margine/exercise1/epw/"
    dir = "./run/epw/"

    wsph, a2f_iso = read_a2f(prefix; dir=dir)

    Ec = 0.1

    # compute_invR0(0.00044, Ec, wsph, a2f_iso)
    # compute_invR0(0.00042, Ec, wsph, a2f_iso)
    # compute_invR0(0.0004, Ec, wsph, a2f_iso)

    # compute_λ(0.00044, Ec, wsph, a2f_iso)
    # compute_λ(0.00042, Ec, wsph, a2f_iso)
    # compute_λ(0.0004, Ec, wsph, a2f_iso)

    N = 9
    lnbetas = zeros(Float64, N)
    invR0s = zeros(Float64, N)
    lamus = zeros(Float64, N)

    for i in 1:N
        TinK = 4.9 + 0.1 * i
        T = TinK / ev2Kelvin
        lamus[i] = compute_λ(T, Ec, wsph, a2f_iso)
        invR0s[i] = compute_invR0(T, Ec, wsph, a2f_iso)
        lnbetas[i] = log10(1 / TinK)
    end

    println(lnbetas)
    println(invR0s)
    println(lamus)

    b, k = linreg(lnbetas, invR0s)
    println("k=$k, b=$b, Tc=$(10^(b/k))")

    b, k = linreg(lnbetas, lamus .- 1.0)
    println("k=$k, b=$b, Tc=$(10^(b/k))")

end