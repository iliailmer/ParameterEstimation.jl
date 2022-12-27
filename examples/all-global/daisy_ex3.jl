import ParameterEstimation
using ModelingToolkit, DifferentialEquations#, Plots
solver = Tsit5()

@parameters p1 p3 p4 p6 p7
@variables t x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
D = Differential(t)

ic = [1.0, 2.0, 1.0, 1]
time_interval = [0.0, 10.0]
datasize = 20
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.2, 0.3, 0.5, 0.6, -0.2] # True Parameters

states = [x1, x2, x3, u0]
parameters = [p1, p3, p4, p6, p7]
@named model = ODESystem([D(x1) ~ -1 * p1 * x1 + x2 + u0,
                             D(x2) ~ p3 * x1 - p4 * x2 + x3,
                             D(x3) ~ p6 * x1 - p7 * x3,
                             D(u0) ~ 1],
                         t, states, parameters)
measured_quantities = [y1 ~ x1 + x3, y2 ~ x2]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval; solver = solver)
println(res)