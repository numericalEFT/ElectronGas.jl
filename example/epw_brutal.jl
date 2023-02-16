
# use brutal force way to compute Tc from EPW data

include("epw_io.jl")

const ev2Kelvin = 1.160451812e4

fermi_mat(n, β, N) = π * (2(n - (N ÷ 2) - 1) + 1) / β

function lambdar_iso_fake(w, wsph, a2f_iso)
    g = 0.2
    Ω = 0.005
    return g * Ω^2 / (Ω^2 + w^2)
end

function s_matrix(N, wsph, a2f_iso, β;
    muc=0.1, zcorrection=true, lambdar_func=lambdar_iso)
    # println("$N, $(fermi_mat(1, β, N)), $(fermi_mat(N÷2, β, N)), $(fermi_mat(N, β, N)) ")
    smat = zeros(Float64, (N, N))
    for i in 1:N
        for j in 1:N
            wi, wj = fermi_mat(i, β, N), fermi_mat(j, β, N)
            λ = lambdar_func(wi - wj, wsph, a2f_iso)
            smat[i, j] = λ - muc
            if i == j && zcorrection
                for k in 1:N
                    wk = fermi_mat(k, β, N)
                    smat[i, j] = smat[i, j] - lambdar_func(wi - wk, wsph, a2f_iso) * sign(wi * wk)
                end
            end
            smat[i, j] = smat[i, j] / abs(wi) * π / β
        end
    end
    return smat
end

sqnorm(x) = sqrt(sum(x .^ 2))

function power_method(smat; conv=1e-6, Nmax=200, shift=2.0)
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

    # source = zeros(Float64, size(smat)[1])
    # source[N÷2] = 1
    # source[N÷2+1] = 1
    source = [1 / abs(2(i - (N ÷ 2) - 1) + 1) for i in 1:N]

    invR0 = 0.0
    diff = 1.0
    n = 1
    xsum = source
    x = xsum .* (1 - α)

    converge = false

    while (!converge && n < Nmax)
        xsum = xsum .* α .+ (source .+ smat * x)
        x = xsum .* (1 - α)
        # diff = abs(1 / x[N÷2] - invR0)
        converge = isapprox(1 / x[N÷2], invR0, rtol=conv, atol=1e-10)
        invR0 = 1 / x[N÷2]
        n = n + 1
    end

    return invR0
end

function compute_λ(T, Ec, wsph, a2f_iso; printlv=1, kwargs...)
    β = 1 / T
    N = floor(Int, Ec / π / T / 2) * 2
    smat = s_matrix(N, wsph, a2f_iso, β; kwargs...)
    λ = power_method(smat)
    if printlv == 1
        println("λ=", λ, ", at T=", ev2Kelvin / β)
        # println("λ from LA: $(eigmax(smat))")
    end
    return λ
end

function compute_invR0(T, Ec, wsph, a2f_iso; printlv=1, kwargs...)
    β = 1 / T
    N = floor(Int, Ec / π / T / 2) * 2
    smat = s_matrix(N, wsph, a2f_iso, β; kwargs...)
    invR0 = precursory_cooper_flow(smat)
    if printlv == 1
        println("1/R0=", invR0, ", at T=", ev2Kelvin / β)
    end
    return invR0
end

function linreg(X, Y)
    hcat(fill!(similar(X), 1), X) \ Y
end

function search_tc_pm(Ec, wsph, a2f_iso; T1=0.1Ec, T2=0.001Ec, Nmax=20, kwargs...)
    lam1 = compute_λ(T1, Ec, wsph, a2f_iso; printlv=0, kwargs...)
    lam2 = compute_λ(T2, Ec, wsph, a2f_iso; printlv=0, kwargs...)

    @assert lam1 < 1.0 # not SC at T1
    n = 1
    while lam2 < 1.0 && n < Nmax
        # ensure SC at T2
        T2 = T2 / 10
        lam2 = compute_λ(T2, Ec, wsph, a2f_iso; printlv=0, kwargs...)
        n = n + 1
    end
    if lam2 < 1.0
        println("No SC!")
        return -1.0
    end

    for i in 1:Nmax
        # search Tc
        T3 = sqrt(T1 * T2)
        lam3 = compute_λ(T3, Ec, wsph, a2f_iso; printlv=0, kwargs...)
        if lam3 < 1.0
            lam1 = lam3
            T1 = T3
        else
            lam2 = lam3
            T2 = T3
        end
    end

    return sqrt(T2 * T1)
end

function search_tc_pcf(Ec, wsph, a2f_iso;
    T1=0.05Ec, T2=0.001Ec, Nmax=20, atol=1e-7, kwargs...)
    lam1 = compute_invR0(T1, Ec, wsph, a2f_iso; printlv=0, kwargs...)
    lam2 = compute_invR0(T2, Ec, wsph, a2f_iso; printlv=0, kwargs...)

    @assert abs(lam1) > atol # not SC at T1
    n = 1
    while abs(lam2) > atol && n < Nmax
        # ensure SC at T2
        T2 = T2 / 10
        lam2 = compute_invR0(T2, Ec, wsph, a2f_iso; printlv=0, kwargs...)
        n = n + 1
    end
    if abs(lam2) > atol
        println("No SC!")
        return -1.0
    end

    for i in 1:Nmax
        # search Tc
        T3 = sqrt(T1 * T2)
        lam3 = compute_invR0(T3, Ec, wsph, a2f_iso; printlv=0, kwargs...)
        if abs(lam3) > atol
            lam1 = lam3
            T1 = T3
        else
            lam2 = lam3
            T2 = T3
        end
    end

    return sqrt(T2 * T1)
end

reflectkwargs(; kwargs...) = kwargs

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

    N = 15
    lnbetas = zeros(Float64, N)
    invR0s = zeros(Float64, N)
    lamus = zeros(Float64, N)

    kwargs = reflectkwargs(muc=0.0, zcorrection=false, lambdar_func=lambdar_iso_fake)
    println(kwargs)

    for i in 1:N
        # TinK = 5.0 + 0.2 * (i - 1)
        TinK = 0.2 * sqrt(2)^(i - 1)
        T = TinK / ev2Kelvin
        lamus[i] = compute_λ(T, Ec, wsph, a2f_iso; kwargs...)
        invR0s[i] = compute_invR0(T, Ec, wsph, a2f_iso; kwargs...)
        lnbetas[i] = log10(1 / TinK)
    end

    println(lnbetas)
    println(invR0s)
    println(lamus)
    isplot = true
    if isplot
        using Plots
        p = plot(lnbetas, invR0s, label="1/R0")
        plot!(p, lnbetas, lamus .- 1.0, label="λ-1")
        plot!(p, lnbetas, lnbetas .* 0.0)
        display(p)
        readline()
    end

    b, k = linreg(lnbetas, invR0s)
    println("k=$k, b=$b, Tc=$(10^(b/k))")
    # Tcfind = search_tc_pcf(Ec, wsph, a2f_iso; kwargs...) * ev2Kelvin
    # println("Tc_find=$(Tcfind)")

    b, k = linreg(lnbetas, lamus .- 1.0)
    println("k=$k, b=$b, Tc=$(10^(b/k))")
    # Tcfind = search_tc_pm(Ec, wsph, a2f_iso; kwargs...) * ev2Kelvin
    # println("Tc_find=$(Tcfind)")

end