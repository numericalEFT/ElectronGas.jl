@testset "Parameter" begin
	# test constructor and reconstruct of Para

	beta = 1e4
	rs = 2.0
	param = Parameter.defaultUnit(1 / beta, rs)

	@test param.me == 0.5
	@test param.EF == 1.0
	@test param.kF == 1.0
	@test param.β == beta
	@test param.gs == 1.0
	@test param.ga == 0.0
	@test param.μ ≈ 1.0

	@test param.ωp ≈ sqrt(4 * param.kF^3 * param.e0^2 / 3 / param.me / π)
	@test param.qTF ≈ sqrt(4 * π * param.e0^2 * param.NF)

	newbeta = 1e3
	ga = 1.0
	newparam = Parameter.derive(param, β = newbeta, ga = ga)

	@test newparam.me == 0.5
	@test newparam.EF == 1.0
	@test newparam.kF == 1.0
	@test newparam.β == newbeta
	@test newparam.ga == ga
	@test newparam.μ ≈ 1.0

	@test newparam.ωp ≈ sqrt(4 * newparam.kF^3 * newparam.e0^2 / 3 / newparam.me / π)
	@test newparam.qTF ≈ sqrt(4 * π * newparam.e0^2 * newparam.NF)

	beta = 2.0
	param = Parameter.rydbergUnit(1 / beta, rs)
	@test param.beta == 2.0
	@test param.kF == (9π / 4)^(1 / 3) / rs
	@test param.EF == param.kF^2
	@test param.μ ≈ param.EF * 0.743112084259
end
