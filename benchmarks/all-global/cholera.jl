using ParameterEstimation
using ModelingToolkit, DifferentialEquations

solver = AutoTsit5(Rosenbrock23())

@parameters mu bi bw al g dz k
@variables t s(t) i(t) w(t) r(t) y1(t) y2(t)
D = Differential(t)

states = [s, i, w, r]
parameters = [mu, bi, bw, al, g, dz, k]

@named model = ODESystem([
                             D(s) ~ mu - bi * s * i - bw * s * w - mu * s + al * r,
                             D(i) ~ bw * s * w + bi * s * i - g * i - mu * i,
                             D(w) ~ dz * (i - w),
                             D(r) ~ g * i - mu * r - al * r,
                         ], t, states, parameters)
measured_quantities = [y1 ~ i, y2 ~ i + r + s]

ic = [1.0, -1.0, 1.0, -1.0]
time_interval = [0.0, 1.0]
datasize = 10
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5] # True Parameters
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic,
                                              datasize; solver = solver)
# plot(data_sample[i], label = "data")
# plot!(data_sample[i + r + s], label = "data")

identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
interpolation_degree = 7
res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)
filtered = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                data_sample, time_interval)

res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval)
print(res)
