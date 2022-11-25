import ParameterEstimation

using ModelingToolkit, Plots, Nemo, HomotopyContinuation

@parameters g a b
@variables t V(t) R(t) y1(t)
D = Differential(t)

u0 = [1.0, -1.0]
time_interval = (0.0, 1.2)
datasize = 16
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [2, 2 / 10, 2 / 10] # True Parameters
measured_quantities = [y1 ~ V]
states = [V, R]
parameters = [g, a, b]

@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ])

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