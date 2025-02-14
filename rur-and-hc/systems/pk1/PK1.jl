# PK1
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L542C16-L549C10

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters k1 k2 k3 k4 k5 k6 k7 s2 s3 u1
@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3, x4]
parameters = [k1, k2, k3, k4, k5, k6, k7, s2, s3, u1]

@info "PK1: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(x1) ~ u1 - (k1 + k2) * x1,
            D(x2) ~ k1 * x1 - (k3 + k6 + k7) * x2 + k5 * x4,
            D(x3) ~ k2 * x1 + k3 * x2 - k4 * x3,
            D(x4) ~ k6 * x2 - k5 * x4,
         ], t, states, parameters)
measured_quantities = [y1 ~ s2 * x2, y2 ~ s3 * x3]

ic = [1.0 for i in 1:length(states)]
time_interval = [0.0, 1.0]
datasize = 21
p_true = [0.1 + i*(1 / (2length(parameters))) for i in 1:length(parameters)]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)

@info "" data_sample
println(parameters)
println(p_true)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
        solver = solver)

