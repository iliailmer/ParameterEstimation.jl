using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()
@parameters k1 k2 k3 k4 k5 k6
@variables t x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3, x4, x5, x6]
parameters = [k1, k2, k3, k4, k5, k6]

@named model = ODESystem([
                             D(x1) ~ -k1 * x1 * x2 + k2 * x4 + k4 * x6,
                             D(x2) ~ -k1 * x1 * x2 + k2 * x4 + k3 * x4,
                             D(x3) ~ k3 * x4 + k5 * x6 - k6 * x3 * x5,
                             D(x4) ~ k1 * x1 * x2 - k2 * x4 - k3 * x4,
                             D(x5) ~ k4 * x6 + k5 * x6 - k6 * x3 * x5,
                             D(x6) ~ -k4 * x6 - k5 * x6 + k6 * x3 * x5,
                         ], t, states, parameters)
measured_quantities = [y1 ~ x3, y2 ~ x2]

ic = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
p_true = [0.03, 0.02, 0.05, 0.03, 0.02, 0.05] # True Parameters

time_interval = [0.0, 30.0]
datasize = 10

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)

res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval; solver = solver)
println(res)