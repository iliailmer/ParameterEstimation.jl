using ParameterEstimation
using ModelingToolkit # ODE definitions

# define toy model
@parameters mu
@variables t x1(t) y1(t) u(t) [input = true]
D = Differential(t)
@named model = ODESystem([D(x1) ~ -mu * x1],
                         t, [x1], [mu])
# data sampling
outputs = [y1 ~ x1 + x1^2]
time = [0.0, 1] # data sample interval

data_sample = Dict{Any, Vector{Float64}}(x1^2 + x1 => [2.0, 1.98506, 1.17611, 0.97441],
                                         "t" => [0.0, 0.01, 0.73, 1.0])

# interpolation and estimation (steps 1-5)
res = ParameterEstimation.estimate(model, outputs, data_sample, time)
