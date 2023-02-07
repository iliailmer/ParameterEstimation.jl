using ParameterEstimation
using ModelingToolkit, DifferentialEquations#, Plots
solver = Tsit5()

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = [0.0, 1.0]
datasize = 20
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
                         states, parameters)

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
# ParameterEstimation.write_sample(data_sample;
#  filename = "../matlab/amigo_models/lotka-volterra-$datasize-$(time_interval[1])-$(time_interval[end]).txt")

# res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
#                                    solver = solver)
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
        title = "Lotka-Volterra Model, $(num_unknowns) unknowns", legend = false)
plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
      title = "Lotka-Volterra Model, $(num_unknowns) unknowns", legend = false)

png("lv_t_$(time_interval[1])_$(time_interval[2]).png")