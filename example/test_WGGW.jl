using LinearAlgebra
using MCIntegration

const T = 0.001
const nmax = 200000
const kmax = 10.0

ε(k) = k^2 / 2 - 1.0
Π(k, ω) = 1 / (ε(k)^2 + ω^2)
Ws(q, ω; Λ=0.0, ωp=1.0) = 4π / (q^2 + Λ) * ω^2 / (ω^2 + ωp^2)
# Wsavg(q, ω; Λ=0.0, ωp=1.0) = 4π * log(q^2 + Λ) * ω^2 / (ω^2 + ωp^2)

# integrand(k, n; T=T) = T * k^2 * Wsavg(k - 1, π * T * (2n))^2 * Π(k, π * T * (2n + 1))

function integrand2(k, n, θ; T=T)
    # integrate over k2, θ, n
    # T\sum_n\int dk2^3/(2π)^3 \int_θ sin(θ) Ws(k1-k2, wn)GG(k2, wn+πT)Ws(k2-k3, wn) 
    k1 = [1.0, 0.0, 0.0]
    k2 = [k[1], k[2], k[3]]
    k3 = [cos(θ), sin(θ), 0.0]
    q1 = dot(k1 - k2, k1 - k2)
    q2 = dot(k2 - k3, k2 - k3)
    ksq = dot(k2, k2)
    factor = 1 / (2π)^3
    return factor * T * sin(θ) * Ws(sqrt(q1), π * T * (2n)) * Π(sqrt(ksq), π * T * (2n + 1)) * Ws(sqrt(q2), π * T * (2n))
end

# f((n, x), c) = integrand(x[1], n[1])

# result = integrate(f; var=(Discrete(-nmax, nmax), Continuous(0.0, kmax; alpha=3.0, adapt=true)), dof=[[1, 1],], neval=1e7)
# println(result)

f2((n, x, θ), c) = integrand2(x, n[1], θ[1])
result = integrate(f2; var=(Discrete(100, nmax), Continuous(0.0, kmax; alpha=3.0, adapt=true), Continuous(-π, π; alpha=1.0, adapt=true)), dof=[[1, 3, 1],], neval=1e6)
println(result)
