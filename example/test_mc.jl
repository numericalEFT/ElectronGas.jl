using Parameters
using Plots
using Random

@with_kw struct Para
    g::Float64 = 0.7 / π
    T::Float64 = 2e-2
    E_c::Float64 = π
    Ω::Float64 = 0.1
    f::Float64 = 1.0
end

para = Para()
print(para)

function Veffsim(ω1, ω2; para=para)
    ωsqp = (ω1 + ω2)^2
    ωsqm = (ω1 - ω2)^2
    # return para.g * π * (ωsqp / (ωsqp + para.Ω^2) + ωsqm / (ωsqm + para.Ω^2))
    return para.g * π * (2*para.f - para.Ω^2*(1 / (ωsqp + para.Ω^2) + 1 / (ωsqm + para.Ω^2)))
end

println("Veffsim(1.0, 2.0)=$(Veffsim(1.0, 2.0))")

function integrand(n, m, para, Δm)
    ωn, ωm = π * para.T * (2n + 1), π * para.T * (2m + 1)
    if abs(ωm) < para.E_c
        return -para.T * Veffsim(ωn, ωm; para=para) * Δm / abs(ωm)
    else
        return 0.0
    end
end

function calcΔ_num(Δ; para=para)
    delta = similar(Δ)

    for n in 1:length(delta)
        delta[n] = 0.0
        for m in 1:length(Δ)
            delta[n] += integrand(n, m, para, Δ[m])
        end
    end
    delta .+= 1.0
    return delta
end

# println(calcΔ_num(ones(100)))

function iterate_num(; para=para, Niter=1000, rtol=1e-6, α=0.8)
    n = floor(Int, para.E_c / para.T)
    Δ = ones(n)

    λ, λold = 0.0, 0.5
    for i in 1:Niter
        delta = calcΔ_num(Δ; para=para)
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

        curr=State([1,1], 1)
        count=n
        reweight=ones(2)
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
    it.curr=prop
end

function reweight!(it)
    it.reweight[2] = abs(sum(it.Δhist) / sum(it.Δ0hist))
end

function Δ(it::MC_Iterator, m)
    return (it.Δ0hist[m]*it.reweight[1] + it.Δhist[m]*it.reweight[2])/it.count*length(it.Δhist)
end

function eval(it::MC_Iterator, st::State=it.curr)
    if st.order==1
        return 1/it.reweight[1] / length(it.Δhist)
    else
        return 1/it.reweight[2] * integrand(st.var[1], st.var[2], it.para, Δ(it, st.var[2]))
    end
end

function propose(it::MC_Iterator)
    curr = it.curr

    if rand()<0.2
        # change order
        prop = State(curr.var, (curr.order%2)+1)
    else
        newn = floor(Int, rand()*length(it.Δhist))+1
        if rand()<0.5
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
    return abs(eval(it, prop)/eval(it))
end

function iterate!(it::MC_Iterator; Niter=1e6, reweightn=0)
    for i in 1:Niter
        prop = propose(it)
        if rand()<acc_ratio(it, prop)
            update!(it, prop)
        else
            update!(it)
        end
        if reweightn != 0 && i%reweightn == 0
            reweight!(it)
        end

    end
end

function flow(lnbetas; para=para, iter=iterate_num)
    lams = zeros(length(lnbetas))
    for i in 1:length(lnbetas)
        beta = 10^(lnbetas[i])
        param = Para(T=1 / beta)
        Δ, λ = iter(para=param)
        println("λ=$(λ)")
        lams[i] = λ
    end
    x, y = [lnbetas'; lnbetas' .* 0.0 .+ 1.0]', lams
    b = (x'x) \ (x'y)
    println(b)
    Tc = 10^(-(1 - b[2]) / b[1])
    println("Tc=$(Tc)")

    return lnbetas, lams, b, Tc
end

lnbetas = [2.6:0.1:3.2;]
plt = plot(xlims=(0.0, 4.0))
for i in 1:length(lnbetas)
    beta = 10^(lnbetas[i])
    # critical at 3.0<f<3.5
    param = Para(T=1 / beta, f=3.2)
    delta, λ = iterate_num(para=param)
    println("λ=$(λ)")
    ωn = [π*param.T*(2*n-1) for n in 1:length(delta)]
    plot!(plt, ωn, delta)
    # display(plt)
    # readline()
end
display(plt)
readline()

# delta, λ = iterate_num(para=para)
# println("λ=$λ")
# println(delta)


# it=MC_Iterator(para=para)
# iterate!(it; Niter = 1e7, reweightn = 1e3)
# λ = 1-1/Δ(it, 1)
# println("λ=$λ")
# println(it.Δ0hist)
# println(it.Δhist)


# lnbetas = [2.6:0.1:3.4;]
# lnbetas, lams, b, Tc = flow(lnbetas)
# plt = plot(lnbetas, lams)
# display(plt)
# readline()



