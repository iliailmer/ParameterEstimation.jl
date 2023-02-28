using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters g a b
@variables t V(t) R(t) y1(t) y2(t)
D = Differential(t)
states = [V, R]
parameters = [g, a, b]

ic = [1.0, -1.0]
time_interval = [0.0, 1.0]
datasize = 50
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [2, 2 / 10, 2 / 10] # True Parameters
measured_quantities = [y1 ~ V]

@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ], t, states, parameters)

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
                                   solver = solver)
all_params = vcat(ic, p_true)
for each in res
    estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
    println("Max abs rel. err: ",
            maximum(abs.((result_ode.u .- all_params) ./ (all_params))))
end
num_unknowns = length(ic) + length(p_true)
all_params = vcat(ic, p_true)
using OrderedCollections
size_err_map = OrderedDict{Int, Float64}()
for datasize in 3:31
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                                  p_true, ic, datasize; solver = solver)
    # res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
    #                                    solver = solver)
    ParameterEstimation.write_sample(data_sample;
                                     filename = "point_error_data/samples/fhn/fhn-$datasize.txt")

    # for each in res
    #     estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
    #     size_err_map[datasize] = maximum(100 *
    #                                      abs.((estimates .- all_params) ./ (all_params)))
    # end
end

# using Plots
# scatter(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
#         title = "Fitzhugh-Nagumo Model, $(num_unknowns) unknowns", legend = false)
# plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
#       title = "Fitzhugh-Nagumo Model, $(num_unknowns) unknowns", legend = false)

# png("point_error_data/figures/fitzHughNagumo_t_$(time_interval[1])_$(time_interval[2]).png")
# open("point_error_data/data/fitzHughNagumo_t_$(time_interval[1])_$(time_interval[2]).txt",
#      "w") do f
#     for (k, v) in size_err_map
#         println(f, k, " ", v)
#     end
# end
