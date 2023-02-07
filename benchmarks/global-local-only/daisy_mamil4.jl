import ParameterEstimation
using ModelingToolkit, DifferentialEquations#, Plots
solver = Tsit5()

@parameters k01 k12 k13 k14 k21 k31 k41
@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)

ic = [1.0, 2.0, 1.0, -1.0]
time_interval = [0.0, 1.0]
datasize = 10
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.2, 0.3, 0.5, 0.6, 0.2, 1.1, 0.2] # True Parameters

states = [x1, x2, x3, x4]
parameters = [k01, k12, k13, k14, k21, k31, k41]
@named model = ODESystem([
                             D(x1) ~ -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 -
                                     k31 * x1 - k41 * x1,
                             D(x2) ~ -k12 * x2 + k21 * x1,
                             D(x3) ~ -k13 * x3 + k31 * x1,
                             D(x4) ~ -k14 * x4 + k41 * x1],
                         t, states, parameters)
measured_quantities = [y1 ~ x1, y2 ~ x2, y3 ~ x3 + x4]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
# ParameterEstimation.write_sample(data_sample;
#                                  filename = "./benchmarks/matlab/amigo_models/daisy_mamil4_loc-$datasize.txt")
# res = ParameterEstimation.estimate(model, measured_quantities, data_sample)
# all_params = vcat(ic, p_true)
# for each in res
#     estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
#     println("Max abs rel. err: ", maximum(abs.((estimates - all_params) ./ all_params)))
# end

num_unknowns = length(ic) + length(p_true)
all_params = vcat(ic, p_true)
using OrderedCollections
size_err_map = OrderedDict{Int, Float64}()
for datasize in 3:21
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                                  p_true, ic, datasize; solver = solver)

    res = ParameterEstimation.estimate(model, measured_quantities, Dict(data_sample);
                                       solver = solver)
    ParameterEstimation.write_sample(data_sample;
                                     filename = "point_error_data/samples/daisy_mamil4/daisy_mamil4-$datasize.txt")
    for each in res
        estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
        size_err_map[datasize] = maximum(100 *
                                         abs.((estimates .- all_params) ./ (all_params)))
    end
end

using Plots
scatter(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
        title = "DAISY 4-Compartment Model, $(num_unknowns) unknowns", legend = false)
plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
      title = "DAISY 4-Compartment Model, $(num_unknowns) unknowns", legend = false)

png("point_error_data/figures/daisy_mamil4_t_$(time_interval[1])_$(time_interval[2]).png")
open("point_error_data/data/daisy_mamil4_$(time_interval[1])_$(time_interval[2]).txt",
     "w") do f
    for (k, v) in size_err_map
        println(f, "$k $v")
    end
end
