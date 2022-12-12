using ModelingToolkit, DifferentialEquations, Plots

using ParameterEstimation

@parameters a b c d
@variables t x1(t) x2(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2]
parameters = [a, b, c, d]

@named model = ODESystem([
                             D(x1) ~ a * x1 + b * x2,
                             D(x2) ~ c * x1 + d * x2,
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ x2,
]

ic = [1.0, -1.0]
time_interval = (0.0, 10.0)
datasize = 10
tdata_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                               p_true, u0,
                                               datasize; solver = solver)
identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
interpolation_degree = 7
res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)
filtered = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                data_sample, time_interval)

res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval)
