using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters k5 k6 k7 k8 k9 k10
@variables t x4(t) x5(t) x6(t) x7(t) y1(t) y2(t)
D = Differential(t)
states = [x4, x5, x6] #, x7]
parameters = [k5, k6, k7, k8, k9, k10]

@named model = ODESystem([
                             D(x4) ~ -k5 * x4 / (k6 + x4),
                             D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6),
                             D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
                             #  D(x7) ~ k9 * x6 * (k10 - x6) / k10,
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x4,
    y2 ~ x5,
]

ic = [1.0, 1.0, 1.0]
time_interval = [0.0, 1.0]
datasize = 20
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1] # True Parameters

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
ParameterEstimation.write_sample(data_sample;
                                 filename = "./matlab/amigo_models/biohydrogenation-$datasize-$(time_interval[1])-$(time_interval[2]).txt")

res = ParameterEstimation.estimate(model, measured_quantities, data_sample,
                                   time_interval)
print(res)