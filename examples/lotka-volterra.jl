import ParameterEstimation

using ModelingToolkit, DifferentialEquations, Plots
using Nemo, HomotopyContinuation
using TaylorSeries

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

u0 = [100.0, 100.0]
time_interval = (0.0, 1.0)
datasize = 15
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.03, 0.02, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ -k3 * w + k2 * r * w])

prob_true = ODEProblem(model, u0, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, Tsit5(), p = p_true, saveat = tsteps)
data_sample = solution_true[1, :]

interpolation_degree = 7
results = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                       time_interval,
                                       interpolation_degree)
print(results)
best_result = ParameterEstimation.filter_solutions(results, model, time_interval,
                                                   data_sample)