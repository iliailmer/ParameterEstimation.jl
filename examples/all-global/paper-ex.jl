using ParameterEstimation
using ModelingToolkit

# Input:
# -- Differential model
@parameters mu
@variables t x(t) y(t)
D = Differential(t)
@named Σ = ODESystem([D(x) ~ -mu * x],
                     t, [x], [mu])
outs = [y ~ x + x^2]

# -- Data
time = [0.0, 1.0] # sampling interval
data = Dict("t" => [0.000, 0.333, 0.666, 1.000],
            x + x^2 => [2.000, 1.563, 1.229, 0.974])

# Run
res = estimate(Σ, outs, data)

identifiability_result = ParameterEstimation.check_identifiability(Σ;
                                                                   measured_quantities = outs)
interpolation_degree = 2
res = ParameterEstimation.estimate_fixed_degree(Σ, outs, data,
                                                identifiability_result,
                                                interpolation_degree)
ParameterEstimation.filter_solutions(res, identifiability_result, Σ, data)