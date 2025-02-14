# NFkB
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd40.158f3e0.1f816bf0.149cf0.1546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L1081

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters i1 i1a k_prod t1 t2
@variables t u(t) x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t) x11(t) x12(t) x13(t) x14(t) x15(t) y1(t) y2(t) y3(t) y4(t) y5(t) y6(t) y7(t)
D = Differential(t)
states = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, u]
parameters = [i1, i1a, k_prod, t1, t2]

@info "NFkB: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(x1) ~ k_prod - 0.1 * x1 - 0.1 * x1 * u,
            D(x2) ~
                -0.1 * x2 - 0.1 * x2 - 0.1 * x2 * x10 + t1 * x4 -
                0.1 * x2 * x13 +
                t2 * x5 +
                (0.1 * x1 - 0.1 * x2 * x8) * u,
            D(x3) ~ 0.1 * x2 - 0.1 * x3 + 0.1 * x2 * x8 * u,
            D(x4) ~ 0.1 * x2 * x10 - t1 * x4,
            D(x5) ~ 0.1 * x2 * x13 - t2 * x5,
            D(x6) ~ 0.1 * x13 - 0.1 * x6 * x10 + t2 * x5 - i1 * x6,
            D(x7) ~ i1 * 0.1 * x6 - 0.1 * x11 * x7,
            D(x8) ~ 0.1 * x9 - 0.1 * x8,
            D(x9) ~ 0.1 - 0.1 * x7 - 0.1 * x9,
            D(x10) ~
                -0.1 * x2 * x10 - 0.1 * x10 * x6 + 0.1 * x12 - 0.1 * x10 -
                i1a * x10 + 0.1 * x11,
            D(x11) ~ -0.1 * x11 * x7 + i1a * 0.1 * x10 - 0.1 * 0.1 * x11,
            D(x12) ~ 0.1 + 0.1 * x7 - 0.1 * x12,
            D(x13) ~ 0.1 * x10 * x6 - 0.1 * x13 - 0.1 * x2 * x13 + 0.1 * x14,
            D(x14) ~ 0.1 * x11 * x7 - 0.1 * 0.1 * x14,
            D(x15) ~ 0.1 + 0.1 * x7 - 0.1 * x15,
            D(u) ~ 1
         ], t, states, parameters)
measured_quantities = [
    y1 ~ x7,
    y2 ~ x10 + x13,
    y3 ~ x9,
    y4 ~ x1 + x2 + x3,
    y5 ~ x2,
    y6 ~ x12,
    y7 ~ u
]

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

