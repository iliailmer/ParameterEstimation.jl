import ParameterEstimation

using ModelingToolkit, DifferentialEquations#, Plots
using Nemo, HomotopyContinuation
solver = AutoTsit5(Rosenbrock23())

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = (0.0, 1.0)
datasize = 15
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.03, 0.02, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
                         states, parameters)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, solver, p = p_true, saveat = tsteps)
data_sample = Dict(r => solution_true[r])

interpolation_degree = 7
identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)
best_result = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                   data_sample, time_interval;
                                                   solver = solver)
res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval)
println(best_result)
println(res)