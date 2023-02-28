using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters k1 k2 k3 k4 k5 k6 Ki
@variables t x1(t) x2(t) x3(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3]
parameters = [k1, k2, k3, k4, k5, k6, Ki]

@named model = ODESystem([
                             D(x1) ~ k1 * Ki^10 / (Ki^10 + x3^10) - k2 * x1,
                             D(x2) ~ k3 * x1 - k4 * x2,
                             D(x3) ~ k5 * x2 - k6 * x3],
                         t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ x3,
]

ic = [1 / 10, 2 / 10, 25 / 10]
time_interval = [0.0, 1.0]
datasize = 50
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1 / 10, 1, 1 / 10, 1, 1 / 10, 1] # True Parameters
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
                                   method = :msolve)
