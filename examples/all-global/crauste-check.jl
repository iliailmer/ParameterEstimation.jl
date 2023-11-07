using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = AutoTsit5(Rosenbrock23())

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
		D(N) ~ -N * mu_N - N * P * delta_NE,
		D(E) ~ N * P * delta_NE - E^2 * mu_EE -
			   E * delta_EL + E * P * rho_E,
		D(S) ~ S * delta_EL - S * delta_LM - S^2 * mu_LL -
			   E * S * mu_LE,
		D(M) ~ S * delta_LM - mu_M * M,
		D(P) ~ P^2 * rho_P - P * mu_P - E * P * mu_PE -
			   S * P * mu_PL,
	], t, states, parameters)
measured_quantities = [y1 ~ N, y2 ~ E, y3 ~ S + M, y4 ~ P]

ic = [1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [-0.5, 0.5]
datasize = 20
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5, 1.0, 1.0, 1.0, 1.0, 0.9, 1.2] # True Parameters
p_alternate = [1.00000000756707, 1.2999995953073495,
	21.345261243944176, 14.168699723072145,
	1.1958082932770018, 0.6315102053491507,
	1.2723560164942402, 1.5119871588063991,
	0.9999999861891795, 1.0000000180382096,
	-0.47276952581306436, 0.8999995611817182, -1.4731385228137435]
ic_alternate = [1.0001253186524368, 0.9997329711131014, -1.3590750984485447, 3.374897457869863, 1.001956206823529]

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver)

data_sample2 = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_alternate, ic_alternate, datasize; solver = solver)

println(data_sample)
println(data_sample2)

#res = ParameterEstimation.estimate(model, measured_quantities, data_sample; solver = solver)
