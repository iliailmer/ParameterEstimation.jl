using ModelingToolkit, DifferentialEquations
using ParameterEstimation
solver = Tsit5()

@parameters lm d beta a k u c q b h
@variables t x(t) y(t) v(t) w(t) z(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [x, y, v, w, z]
parameters = [lm, d, beta, a, k, u, c, q, b, h]

@named model = ODESystem([
                             D(x) ~ lm - d * x - beta * x * v,
                             D(y) ~ beta * x * v - a * y,
                             D(v) ~ k * y - u * v,
                             D(w) ~ c * x * y * w - c * q * y * w - b * w,
                             D(z) ~ c * q * y * w - h * z,
                         ], t, states, parameters)
measured_quantities = [y1 ~ w, y2 ~ z, y3 ~ x, y4 ~ y + v]

ic = [1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [0.0, 10.0]
datasize = 20
p_true = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic,
                                              datasize; solver = solver)

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
        title = "HIV Model, $(num_unknowns) unknowns", legend = false)
plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
      title = "HIV Model, $(num_unknowns) unknowns", legend = false)

png("point_error_data/figures/hiv_t_$(time_interval[1])_$(time_interval[2]).png")
open("point_error_data/data/hiv_t_$(time_interval[1])_$(time_interval[2]).txt", "w") do f
    for (k, v) in size_err_map
        println(f, k, " ", v)
    end
end
