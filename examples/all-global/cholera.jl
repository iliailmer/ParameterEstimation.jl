
using ModelingToolkit, DifferentialEquations, Plots
using Nemo, HomotopyContinuation

import ParameterEstimation

@parameters mu bi bw al g dz k
@variables t s(t) i(t) w(t) r(t) y1(t) y2(t)
D = Differential(t)
@named model = ODESystem([
                             D(s) ~ mu - bi * s * i - bw * s * w - mu * s + al * r,
                             D(i) ~ bw * s * w + bi * s * i - g * i - mu * i,
                             D(w) ~ dz * (i - w),
                             D(r) ~ g * i - mu * r - al * r,
                         ])
measured_quantities = [y1 ~ i, y2 ~ i + r + s]

u0 = [1.0, -1.0, 1.0, -1.0]
time_interval = (0.0, 1.0)
datasize = 50
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5] # True Parameters
states = [s, i, w, r]
parameters = [mu, bi, bw, al, g, dz, k]

prob_true = ODEProblem(model, u0, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true,
                                      ARKODE(Sundials.Explicit(),
                                             etable = Sundials.FEHLBERG_6_4_5), p = p_true,
                                      saveat = tsteps)

data_sample = Dict(Num(v.rhs) => solution_true[Num(v.rhs)] for v in measured_quantities)
plot(solution_true)
identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
interpolation_degree = 10
res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)

res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval)
