using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters mu_N mu_EE mu_LE mu_LL mu_M mu_P mu_PE mu_PL delta_NE delta_EL delta_LM rho_E rho_P
@variables t N(t) E(t) S(t) M(t) P(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [N, E, S, M, P]
parameters = [
	mu_N,
	mu_EE,
	mu_LE,
	mu_LL,
	mu_M,
	mu_P,
	mu_PE,
	mu_PL,
	delta_NE,
	delta_EL,
	delta_LM,
	rho_E,
	rho_P,
]
@named model = ODESystem(
	[
		D(N) ~ -N * mu_N^2 - N * P * delta_NE^2,
		D(E) ~ N * P * delta_NE - E^2 * mu_EE^2 -
			   E * delta_EL + E * P * rho_E^2,
		D(S) ~ S * delta_EL - S * delta_LM - S^2 * mu_LL -
			   E * S * mu_LE,
		D(M) ~ S * delta_LM - mu_M * M,
		D(P) ~ P^2 * rho_P - P * mu_P - E * P * mu_PE -
			   S * P * mu_PL,
	], t, states, parameters)
measured_quantities = [y1 ~ N, y2 ~ E, y3 ~ S + M, y4 ~ P]

ic = [1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [-0.5, 0.5]
datasize = 21
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5, 1.0, 1.0, 1.0, 1.0, 0.9, 1.2] # True Parameters
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)

