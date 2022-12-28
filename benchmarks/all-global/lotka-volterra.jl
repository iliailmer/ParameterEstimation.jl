import ParameterEstimation

using ModelingToolkit, DifferentialEquations#, Plots
solver = Tsit5()

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = [0.0, 1.0]
datasize = 16
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
                         states, parameters)

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
# ParameterEstimation.write_sample(data_sample;
#                                  filename = "benchmarks/matlab/amigo_models/lotka-volterra-$datasize.txt")
interpolation_degree = 9

identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)

best_result = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                   data_sample, time_interval;
                                                   solver = solver)
# print(best_result)
res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval; solver = solver)
println(res)