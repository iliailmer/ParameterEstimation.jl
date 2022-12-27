using ModelingToolkit, DifferentialEquations
using ParameterEstimation
solver = Tsit5()

@parameters mu_N mu_EE mu_LE mu_LL mu_M mu_P mu_PE mu_PL delta_NE delta_EL delta_LM rho_E rho_P
@variables t N(t) E(t) S(t) M(t) P(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [N, E, S, M, P]
parameters = [
    mu_N,
    mu_EE,
    mu_LE,
    mu_LL,
    mu_M,
    mu_P,
    mu_PE,
    mu_PL,
    delta_NE,
    delta_EL,
    delta_LM,
    rho_E,
    rho_P,
]
@named model = ODESystem([
                             D(N) ~ -N * mu_N - N * P * delta_NE,
                             D(E) ~ N * P * delta_NE - E^2 * mu_EE -
                                    E * delta_EL + E * P * rho_E,
                             D(S) ~ S * delta_EL - S * delta_LM - S^2 * mu_LL -
                                    E * S * mu_LE,
                             D(M) ~ S * delta_LM - mu_M * M,
                             D(P) ~ P^2 * rho_P - P * mu_P - E * P * mu_PE -
                                    S * P * mu_PL,
                         ], t, states, parameters)
measured_quantities = [y1 ~ N, y2 ~ E, y3 ~ S + M, y4 ~ P]

ic = [1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [0.0, 10.0]
datasize = 50
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5, 1.0, 1.0, 1.0, 1.0, 0.9, 1.2] # True Parameters
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
# ParameterEstimation.write_sample(data_sample;
#                                  filename = "benchmarks/matlab/amigo_models/crauste-10.txt")

identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = measured_quantities)
interpolation_degree = 8
res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval, identifiability_result,
                                   interpolation_degree)

filtered = ParameterEstimation.filter_solutions(res, identifiability_result, model,
                                                data_sample, time_interval; solver = solver)
print(filtered)
res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval; solver = solver)
print(res)