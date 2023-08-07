using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters a b
@variables t x1(t)
D = Differential(t)
states = [x1, x2]
parameters = [a, b]

@named model = ODESystem([
		D(x1) ~ -a * x1,
		D(x2) ~ -b * x2], t, states, parameters)
measured_quantities = [
	y1 ~ x1,
	y2 ~ x2,
]

ic = [1.0, 1.0]
p_true = [1.0, 1.0]
time_interval = [0.1, 4.0]
datasize = 20
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample)
