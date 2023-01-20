using ParameterEstimation
using ModelingToolkit

# Input:
# -- Differential model
@parameters mu
@variables t x(t) y(t)
D = Differential(t)
@named Sigma = ODESystem([D(x) ~ -mu * x],
                         t, [x], [mu])
outs = [y ~ x + x^2]

# -- Data
time = [0.0, 1.0] # sampling interval
data = Dict("t" => [0.0, 0.333, 0.666, 1.0],
            x + x^2 => [2.0, 1.563, 1.229, 0.974])

# Run
res = estimate(Sigma, outs, data, time)
