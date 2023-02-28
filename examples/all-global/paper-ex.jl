using ParameterEstimation
using ModelingToolkit

# Input:
# -- Differential model
@parameters mu
@variables t x(t) y(t)
D = Differential(t)
<<<<<<< HEAD
@named Sigma = ODESystem([D(x) ~ -mu * x],
                         t, [x], [mu])
outs = [y ~ x + x^2]

# -- Data
||||||| c9f200c
@named model = ODESystem([D(x1) ~ -mu * x1],
                         t, [x1], [mu])
# data sampling
outs = [y1 ~ x1 + x1^2]
=======
@named Σ = ODESystem([D(x) ~ -mu * x],
                     t, [x], [mu])
outs = [y ~ x + x^2]

# -- Data
>>>>>>> main
time = [0.0, 1.0] # sampling interval
<<<<<<< HEAD
data = Dict("t" => [0.0, 0.333, 0.666, 1.0],
            x + x^2 => [2.0, 1.563, 1.229, 0.974])

# Run
res = estimate(Sigma, outs, data, time)
||||||| c9f200c
data = Dict{Any, Vector{Float64}}("t" => [0.0, 1 / 3, 2 / 3, 1.0],
                                  x1 + x1^2 => [2.0, 1.56301,
                                      1.22995, 0.97441])
# identifiabliity,
# interpolation and estimation (steps 1-5)
res = estimate(model, outs, data, time)
identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                   measured_quantities = outs)
res = ParameterEstimation.estimate_fixed_degree(model, outs, data, time,
                                                identifiability_result, 1)
=======
data = Dict("t" => [0.000, 0.333, 0.666, 1.000],
            x + x^2 => [2.000, 1.563, 1.229, 0.974])

# Run
res = estimate(Σ, outs, data)
>>>>>>> main
