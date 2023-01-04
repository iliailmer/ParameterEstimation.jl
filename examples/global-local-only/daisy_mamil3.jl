using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters a21 a31 a01 a12 a13 a31
@variables t x1(t) x2(t) x3(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3]
parameters = [a21, a31, a01, a12, a13, a31]

@named model = ODESystem([
                             D(x1) ~ -(a21 + a31 + a01) * x1 + a12 * x2 +
                                     a13^2 * x3,
                             D(x2) ~ a21 * x1 - a12 * x2,
                             D(x3) ~ a31 * x1 - a13^2 * x3,
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ x2,
]

ic = [1.0, -1.0, 1.0]
time_interval = [0.0, 2.0]
datasize = 21
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1.3, 1.1, 1.2, 1, 1.2] # True Parameters
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)

identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
interpolation_degree = 15
res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)
filtered = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                data_sample, time_interval)
print(filtered)
