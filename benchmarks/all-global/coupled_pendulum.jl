import ParameterEstimation

using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters alpha
@variables t theta_1(t) theta_1_dot(t) theta_2(t) theta_2_dot(t) y1(t) y2(t)
D = Differential(t)

#coupled pendulum equations
states = [theta_1, theta_1_dot, theta_2, theta_2_dot]
parameters = [alpha]
@named model = ODESystem([
                             D(theta_1) ~ theta_1_dot,
                             D(theta_1_dot) ~ -theta_1 * (alpha + 1) + alpha * theta_2,
                             D(theta_2) ~ theta_2_dot,
                             D(theta_2_dot) ~ alpha * theta_1 - theta_2 * (alpha + 1),
                         ], t, states, parameters)

#initial conditions
ic = [1.0, 0.0, 0.0, -1.0]
time_interval = [0.0, 10.0]
datasize = 100
p_true = [1.0]
measured_quantities = [y1 ~ theta_1, y2 ~ theta_2 + theta_1]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
# identifiability_result = ParameterEstimation.check_identifiability(model;
#                                                                    measured_quantities = measured_quantities)
res = ParameterEstimation.estimate_over_degrees(model, measured_quantities, data_sample,
                                                time_interval; solver = solver)
print(res)