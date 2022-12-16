import ParameterEstimation

using ModelingToolkit, Plots, Nemo, HomotopyContinuation, DifferentialEquations

@parameters g a b
@variables t V(t) R(t) y1(t) y2(t)
D = Differential(t)
states = [V, R]
parameters = [g, a, b]

u0 = [1.0, -1.0]
time_interval = (0.0, 5)
datasize = 50
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [2, 2 / 10, 2 / 10] # True Parameters
measured_quantities = [y1 ~ V]

@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ], t, states, parameters)

prob_true = ODEProblem(model, u0, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true,
                                      ARKODE(Sundials.Explicit(),
                                             etable = Sundials.FEHLBERG_6_4_5), p = p_true,
                                      saveat = tsteps)
data_sample = Dict(V => solution_true[V])
at_time = 0.0
# interpolation_degree = 23
# identifiability_result = ParameterEstimation.check_identifiability(model;
#                                                                    measured_quantities = measured_quantities)
# res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
#                                    time_interval, identifiability_result,
#                                    interpolation_degree, at_time)

# filtered = ParameterEstimation.filter_solutions(res, identifiability_result, model,
#                                                 data_sample, time_interval)
res = ParameterEstimation.estimate_over_degrees(model, measured_quantities,
                                                data_sample,
                                                time_interval, at_time)

time_res = Dict()
idx = 1
for t in range(time_interval[1], time_interval[2], length = 10)
    time_res[idx] = ParameterEstimation.estimate_over_degrees(model, measured_quantities,
                                                              data_sample,
                                                              time_interval, t)
    idx += 1
end
