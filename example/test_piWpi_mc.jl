using MCIntegration

Π(ω; ωc=0.1) = 1 / abs(ω) * atan(ωc / abs(ω))
D(ω; ωD=0.001) = ωD^2 / (ω^2 + ωD^2)
Ws(ω; ωp=0.1) = ω^2 / (ω^2 + ωp^2)

f(x, c) = Ws(x[1]) * Π(x[1]) * Π(x[2]) * D(x[1] - x[2]; ωD=0.001) / π^2

T = 0.001
result = integrate(f; var=Continuous(0.882T, 10.0), dof=[[2,],], neval=1e7)
println(result)