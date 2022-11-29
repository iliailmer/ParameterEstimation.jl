# Citation: Zamudio Lara, J.M.; Dewasme, L.; Hernández Escoto, H.; Vande Wouwer, A.
# Parameter Estimation of Dynamic Beer Fermentation Models. Foods 2022, 11, 3602.

using ModelingToolkit, DifferentialEquations, Plots
using Nemo, HomotopyContinuation

import ParameterEstimation

@parameters k5 k6 k7 k8 k9 k10
@variables t x4(t) x5(t) x6(t) x7(t) y1(t) y2(t)
D = Differential(t)
@named model = ODESystem([
                             D(x4) ~ -k5 * x4 / (k6 + x4),
                             D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6),
                             D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
                             #  D(x7) ~ k9 * x6 * (k10 - x6) / k10,
                         ])
measured_quantities = [
    y1 ~ x4,
    y2 ~ x5,
]

u0 = [1.0, -1.0, 1.0]
time_interval = (0.0, 6.0)
datasize = 50
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1] # True Parameters
states = [x4, x5, x6] #, x7]
parameters = [k5, k6, k7, k8, k9, k10]

prob_true = ODEProblem(model, u0, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, Tsit5(), p = p_true, saveat = tsteps)

data_sample = Dict(Num(v.rhs) => solution_true[Num(v.rhs)] for v in measured_quantities)
# plot(solution_true)
identifiability_result = ParameterEstimation.get_identifiability(model;
                                                                 measured_quantities = measured_quantities)
interpolation_degree = 10
results = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                       time_interval, identifiability_result,
                                       interpolation_degree)

results = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                    time_interval)
