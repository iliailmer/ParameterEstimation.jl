
using ModelingToolkit, DifferentialEquations, Plots
using Nemo, HomotopyContinuation

import ParameterEstimation

@parameters p1 p2 p3 p4 p6 p7
@variables t x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
D = Differential(t)
@named model = ODESystem([
                             D(x1) ~ -1 * p1 * x1 + x2 + u0,
                             D(x2) ~ p3 * x1 - p4 * x2 + x3,
                             D(x3) ~ p6 * x1 - p7 * x3,
                             D(u0) ~ 1,
                         ])
measured_quantities = [
    y1 ~ x1,
    y2 ~ u0,
]

ic = [1.0, -1.0, 1.0, -1.0]
time_interval = (0.0, 6.0)
datasize = 50
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1.3, 1.1, 1.2, 1] # True Parameters
states = [x1, x2, x3, u0]
parameters = [p1, p3, p4, p6, p7]

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, Tsit5(), p = p_true, saveat = tsteps)

data_sample = Dict(Num(v.rhs) => solution_true[Num(v.rhs)] for v in measured_quantities)
# plot(solution_true)
identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
interpolation_degree = 37
results = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                       time_interval, identifiability_result,
                                       interpolation_degree)
filtered = ParameterEstimation.filter_solutions(results, identifiability_result, model,
                                                data_sample, time_interval)

results = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                    time_interval)
