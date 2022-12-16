import ParameterEstimation

using ModelingToolkit, DifferentialEquations#, Plots
using Sundials
using Nemo, HomotopyContinuation
solver = Tsit5()

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = (0.0, 1.0)
datasize = 32
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
                         states, parameters)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, solver, p = p_true, saveat = tsteps,
                                      abstol = 1e-12, reltol = 1e-12)

data_sample = Dict(r => solution_true[r])

interpolation_degree = 9
identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
num_parameters = length(parameters) + length(states)
@info "Interpolating sample data"
polynomial_system, interpolants = ParameterEstimation.interpolate(identifiability_result,
                                                                  data_sample,
                                                                  time_interval,
                                                                  measured_quantities,
                                                                  interpolation_degree,
                                                                  num_parameters + 1,
                                                                  0.0)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)

function relative_error(p_true, p_est)
    return abs(p_true - p_est) / abs(p_true)
end
best_result = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                   data_sample, time_interval;
                                                   solver = solver)
maximum(map(x -> relative_error(x...), zip(values(best_result[1].parameters), p_true)))

res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval;
                                                solver = solver)
# # println(best_result)
# println(res)