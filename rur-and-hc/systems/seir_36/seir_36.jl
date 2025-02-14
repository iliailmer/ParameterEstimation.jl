# SEIR_36_ref
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L672

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters beta beta_d gamma gamma_d mu_0 mu_d mu_i phi phi_e s s_d
@variables t N(t) nu(t) q(t) S(t) E(t) I(t) De(t) Di(t) R(t) F(t) y1(t) y2(t) y3(t) y4(t) y5(t) y6(t)
D = Differential(t)
states = [N, nu, q, S, E, I, De, Di, R, F]
parameters = [beta, beta_d, gamma, gamma_d, mu_0, mu_d, mu_i, phi, phi_e, s, s_d]

@info "SEIR_36_ref: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(N) ~ 0,
            D(nu) ~ 0,
            D(q) ~ 0,
            D(S) ~
                -beta * S * I / N - q * beta_d * S * Di / N +
                nu * N - mu_0 * S,
            D(E) ~
                beta * S * I / N + q * beta_d * S * Di / N - s * E - phi_e * E - mu_0 * E,
            D(I) ~ s * E - gamma * I - mu_i * I - phi * I - mu_0 * I,
            D(De) ~ phi_e * E - s_d * De - mu_0 * De,
            D(Di) ~
                phi * I + s_d * De - gamma_d * Di - mu_d * Di - mu_0 * Di,
            D(R) ~ gamma * I + gamma_d * Di - mu_0 * R,
            D(F) ~ mu_i * I + mu_d * Di,
         ], t, states, parameters)
measured_quantities = [
    y1 ~ De,
    y2 ~ Di,
    y5 ~ F,
    y3 ~ N,
    y4 ~ nu,
    y6 ~ q
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


