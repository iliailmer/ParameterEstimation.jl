using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters p1 p3 p4 p6 p7
@variables t x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3, u0]
parameters = [p1, p3, p4, p6, p7]
@named model = ODESystem([
                             D(x1) ~ -1 * p1 * x1 + x2 + u0,
                             D(x2) ~ p3 * x1 - p4 * x2 + x3,
                             D(x3) ~ p6 * x1 - p7 * x3,
                             D(u0) ~ 1,
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ u0,
]

ic = [1.0, -1.0, 1.0, -1.0]
time_interval = [0.0, 1]
datasize = 20
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1.3, 1.1, 1.2, 1] # True Parameters

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
# ParameterEstimation.write_sample(data_sample;
#  filename = "./benchmarks/matlab/amigo_models/daisy_ex3-loc-$datasize-$(time_interval[1])-$(time_interval[2]).txt")

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

    res = ParameterEstimation.estimate(model, measured_quantities, Dict(data_sample);
                                       solver = solver)
    ParameterEstimation.write_sample(data_sample;
                                     filename = "point_error_data/samples/daisy_ex_3/daisy_ex_3-$datasize.txt")
    for each in res
        estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
        size_err_map[datasize] = maximum(100 *
                                         abs.((estimates .- all_params) ./ (all_params)))
    end
end

using Plots
scatter(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
        title = "DAISY 3-Compartment Model 1, $(num_unknowns) unknowns", legend = false)
plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
      title = "DAISY 3-Compartment Model 1, $(num_unknowns) unknowns", legend = false)

png("point_error_data/figures/daisy_ex_3_t_$(time_interval[1])_$(time_interval[2]).png")
open("point_error_data/data/daisy_ex_3_$(time_interval[1])_$(time_interval[2]).txt",
     "w") do f
    for (k, v) in size_err_map
        println(f, "$k $v")
    end
end
