using Parameters
using Plots
using Random

@with_kw struct Para
    g::Float64 = 1.0 # 0.7 / π
    T::Float64 = 2e-2
    E_c::Float64 = π
    Ω::Float64 = π / 10
    f::Float64 = 0.5
end

form = :rs
# form = :chub
para = Para()
print(para)

function Veffsim(ω1, ω2; para=para, para2=para, iscombine=false, form=form)
    if ω1 > para.E_c || ω2 > para.E_c
        return 0.0
    end
    ωsqp = (ω1 + ω2)^2
    ωsqm = (ω1 - ω2)^2
    # return para.g * π * (ωsqp / (ωsqp + para.Ω^2) + ωsqm / (ωsqm + para.Ω^2))

    if form == :rs
        if ω1 < para.Ω && ω2 < para.Ω
            V = 2 * para.g * π * (1 - para.f)
        else
            V = 2 * para.g * π
        end
    elseif form == :chub
        V = para.g * π * (2 * para.f - para.Ω^2 * (1 / (ωsqp + para.Ω^2) + 1 / (ωsqm + para.Ω^2)))

    else
        V = 0.0
    end

    if iscombine == false
        return V
    else
        return V + Veffsim(ω1, ω2; para=para2, iscombine=false, form=form)
    end
end



function integrand(n, m, para, Δm; iscombine=false, para2=para)
    ωn, ωm = π * para.T * (2n - 1), π * para.T * (2m - 1)
    if abs(ωm) < para.E_c
        return -para.T * Veffsim(ωn, ωm; para=para, iscombine=iscombine, para2=para2) * Δm / abs(ωm)
    else
        return 0.0
    end
end

function calcΔ_num(Δ; para=para, iscombine=false, para2=para)
    delta = similar(Δ)

    for n in 1:length(delta)
        delta[n] = 0.0
        for m in 1:length(Δ)
            delta[n] += integrand(n, m, para, Δ[m]; iscombine=iscombine, para2=para2)
        end
    end
    delta .+= 1.0
    return delta
end

# println(calcΔ_num(ones(100)))

function iterate_num(; para=para, Niter=3000, rtol=1e-8, α=0.85, iscombine=false, para2=para)
    n = floor(Int, para.E_c / para.T / π / 2)
    Δ = ones(n)

    λ, λold = 0.0, 0.5
    for i in 1:Niter
        delta = calcΔ_num(Δ; para=para, iscombine=iscombine, para2=para2)
        λ = 1 - 1 / Δ[1]
        diff = abs(λ - λold) / abs(λold)
        if diff < rtol
            break
        else
            λold = λ
            Δ = (Δ * α + delta * (1 - α))
        end
    end

    return Δ, λ
end


struct State
    var::Vector
    order::Int
end

mutable struct MC_Iterator
    para::Para
    count::Int
    Δhist::Vector
    Δ0hist::Vector
    curr::State
    reweight::Vector

    function MC_Iterator(; para=para)
        n = floor(Int, para.E_c / para.T)
        Δhist = zeros(n)
        Δ0hist = ones(n)

        curr = State([1, 1], 1)
        count = n
        reweight = ones(2)
        return new(para, count, Δhist, Δ0hist, curr, reweight)
    end
end

function update!(it::MC_Iterator, prop::State=it.curr)
    if prop.order == 1
        it.count += 1
        it.Δ0hist[prop.var[1]] += 1
    else
        it.Δhist[prop.var[1]] += 1 * sign(eval(it, prop))
    end
    it.curr = prop
end

function reweight!(it)
    it.reweight[2] = abs(sum(it.Δhist) / sum(it.Δ0hist))
end

function Δ(it::MC_Iterator, m)
    return (it.Δ0hist[m] * it.reweight[1] + it.Δhist[m] * it.reweight[2]) / it.count * length(it.Δhist)
end

function eval(it::MC_Iterator, st::State=it.curr)
    if st.order == 1
        return 1 / it.reweight[1] / length(it.Δhist)
    else
        return 1 / it.reweight[2] * integrand(st.var[1], st.var[2], it.para, Δ(it, st.var[2]))
    end
end

function propose(it::MC_Iterator)
    curr = it.curr

    if rand() < 0.2
        # change order
        prop = State(curr.var, (curr.order % 2) + 1)
    else
        newn = floor(Int, rand() * length(it.Δhist)) + 1
        if rand() < 0.5
            # change n
            prop = State([newn, curr.var[2]], curr.order)
        else
            # change m
            prop = State([curr.var[1], newn], curr.order)
        end
    end
    return prop
end

function acc_ratio(it::MC_Iterator, prop::State)
    return abs(eval(it, prop) / eval(it))
end

function iterate!(it::MC_Iterator; Niter=1e6, reweightn=0)
    for i in 1:Niter
        prop = propose(it)
        if rand() < acc_ratio(it, prop)
            update!(it, prop)
        else
            update!(it)
        end
        if reweightn != 0 && i % reweightn == 0
            reweight!(it)
        end

    end
end

function flow(lnbetas; para=para, iter=iterate_num, iscombine=false, para2=para, isplot=false)
    println("calc flow")
    println(para)
    if iscombine
        println(para2)
    end
    if isplot
        plt = plot()
    end
    lams = zeros(length(lnbetas))
    for i in 1:length(lnbetas)
        beta = 10^(lnbetas[i])
        param = Para(para, T=1 / beta)
        param2 = Para(para2, T=1 / beta)
        Δ, λ = iter(para=param; iscombine=iscombine, para2=param2)
        println("λ=$(λ)")
        lams[i] = λ
        if isplot
            ωn = [π * param.T * (2 * n - 1) for n in 1:length(Δ)]
            plot!(plt, ωn, Δ)
        end
    end
    x, y = [lnbetas'; lnbetas' .* 0.0 .+ 1.0]', lams
    b = (x'x) \ (x'y)
    println(b)
    Tc = 10^(-(1 - b[2]) / b[1])
    println("Tc=$(Tc)")

    if isplot
        display(plt)
        readline()
    end
    return lnbetas, lams, b, Tc
end

# lnbetas = [2.6:0.1:3.2;]
# plt = plot(xlims=(0.0, 4.0))
# for i in 1:length(lnbetas)
#     beta = 10^(lnbetas[i])
#     # critical at 3.0<f<3.5
#     param = Para(T=1 / beta, f=1.0)
#     delta, λ = iterate_num(para=param)
#     println("λ=$(λ)")
#     ωn = [π * param.T * (2 * n - 1) for n in 1:length(delta)]
#     plot!(plt, ωn, delta)
#     # display(plt)
#     # readline()
# end
# display(plt)
# readline()

# delta, λ = iterate_num(para=para)
# println("λ=$λ")
# println(delta)


# it=MC_Iterator(para=para)
# iterate!(it; Niter = 1e7, reweightn = 1e3)
# λ = 1-1/Δ(it, 1)
# println("λ=$λ")
# println(it.Δ0hist)
# println(it.Δhist)

function tc_theory(para::Para)
    α = 0.882
    @unpack g, T, E_c, Ω, f = para
    l = log(E_c / Ω)
    L = (1 + g * l) / g / (g * l * f - 1 + f)
    println("L=$L")
    return Ω / α / exp(L)
end

para1 = Para(g=1.0, T=0.1, E_c=1.0, Ω=0.1, f=0.4)
para2 = Para(g=1.0, T=0.1, E_c=1.0 / 2, Ω=0.1 / 2, f=0.4)

println("Veffsim(1.0, 2.0)=$(Veffsim(1.0, 2.0; para=para1, para2=para2, iscombine=true))")
println("Veffsim(0.7, 0.8)=$(Veffsim(0.7, 0.8; para=para1, para2=para2, iscombine=true))")
println("Veffsim(0.2, 0.3)=$(Veffsim(0.2, 0.3; para=para1, para2=para2, iscombine=true))")
println("Veffsim(0.07, 0.08)=$(Veffsim(0.07, 0.08; para=para1, para2=para2, iscombine=true))")
println("Veffsim(0.02, 0.03)=$(Veffsim(0.02, 0.03; para=para1, para2=para2, iscombine=true))")

lnbetas = [2.4:0.01:2.45;]
# lnbetas = [2.4:0.2:3.2;]
lnbetas, lams, b, Tc = flow(lnbetas; para=para1, iscombine=true, para2=para2, isplot=true)
# lnbetas, lams, b, Tc = flow(lnbetas; para=para2, isplot=true)
println("Tc=$Tc, theory=$(tc_theory(para2))")
println("$(b./log(10))")
plt = plot(lnbetas, lams)
display(plt)
readline()



