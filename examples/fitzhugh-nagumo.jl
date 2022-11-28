import ParameterEstimation

using ModelingToolkit, Plots, Nemo, HomotopyContinuation, DifferentialEquations

@parameters g a b
@variables t V(t) R(t) y1(t)
D = Differential(t)

u0 = [1.0, -1.0]
time_interval = (0.0, 1.2)
datasize = 20
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

interpolation_degree = 14
interpolation_point = 0.0
for interpolation_point in range(time_interval[1], time_interval[2], length = 10)
    println("Interpolation point: $interpolation_point")
    results = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                           time_interval, interpolation_degree,
                                           interpolation_point)

    best_estimate, min_error = ParameterEstimation.filter_solutions(results, model,
                                                                    time_interval,
                                                                    data_sample,
                                                                    interpolation_point)

    println("Best estimate: $best_estimate")
    println("Minimum error: $min_error")
end
results = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                       time_interval, interpolation_degree,
                                       interpolation_point)
# results = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
# time_interval)
plot(solution_true)

initial_conditions = [filtered[1][s] for s in ModelingToolkit.states(model)]
parameter_values = [filtered[1][p] for p in ModelingToolkit.parameters(model)]
prob = ModelingToolkit.ODEProblem(model, initial_conditions,
                                  (interpolation_point, time_interval[2]),
                                  parameter_values)
time_step = (time_interval[2] - time_interval[1]) / length(data_sample)
new_length = Int(floor((time_interval[2] - starting_point) / time_step))

ode_solution = ModelingToolkit.solve(prob, Tsit5(), p = parameter_values,
                                     saveat = range(interpolation_point, time_interval[2],
                                                    length = new_length))
plot!(ode_solution)
