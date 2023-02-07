using ModelingToolkit, DifferentialEquations#, Plots
using ParameterEstimation
solver = Tsit5()

@parameters a b
@variables t x1(t) x2(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2]
parameters = [a, b]

@named model = ODESystem([
                             D(x1) ~ -a * x2,
                             D(x2) ~ 1 / b * (x1),
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ x2,
]

ic = [1.0, 1.0]
p_true = [9.8, 1.3]
time_interval = [0.0, 2.0 * pi * sqrt(1.3 / 9.8)]
datasize = 20
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic,
                                              datasize; solver = solver)
# ParameterEstimation.write_sample(data_sample;
#  filename = "../matlab/amigo_models/simple-$datasize-$(time_interval[1])-$(time_interval[end]).txt")
# res = ParameterEstimation.estimate(model, measured_quantities, data_sample)
# all_params = vcat(ic, p_true)
# for each in res
#     estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
#     println("Max abs rel. err: ",
#             maximum(abs.((result_ode.u .- all_params) ./ (all_params))))
# end
num_unknowns = length(ic) + length(p_true)
all_params = vcat(ic, p_true)
using OrderedCollections
size_err_map = OrderedDict{Int, Float64}()
for datasize in 3:21
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                                  p_true, ic, datasize; solver = solver)
    res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
                                       solver = solver)

    for each in res
        estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
        size_err_map[datasize] = maximum(100 *
                                         abs.((estimates .- all_params) ./ (all_params)))
    end
end

using Plots
scatter(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
        title = "Oscillator Model, $(num_unknowns) unknowns", legend = false)
plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
      title = "Oscillator Model, $(num_unknowns) unknowns", legend = false)

png("simple_t_$(time_interval[1])_$(time_interval[2]).png")