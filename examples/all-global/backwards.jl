using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters a b c d
@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [x1]
parameters = [a]

@named model = ODESystem([
		D(x1) ~ a * x1,
		
	], t, states, parameters)
measured_quantities = [
	y1 ~ x1,
	
]

ic = [1.0]
p_true = [2.0]
time_interval = [1.0, 0.0]
datasize = 9
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver)
    println(data_sample)
#res = ParameterEstimation.estimate(model, measured_quantities, data_sample)
