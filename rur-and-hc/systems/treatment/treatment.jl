# Treatment_io
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L516

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters a b d g nu
@variables t S(t) In(t) Tr(t) N(t) y1(t) y2(t)
D = Differential(t)
states = [S, In, Tr, N]
parameters = [a, b, d, g, nu]

@info "Treatment_io: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(S) ~ -b * S * In / N - d * b * S * Tr / N,
            D(In) ~
                b * S * In / N + d * b * S * Tr / N - (a + g) * In,
            D(Tr) ~ g * In - nu * Tr,
            D(N) ~ 0,
         ], t, states, parameters)
measured_quantities = [y1 ~ Tr, y2 ~ N]

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

