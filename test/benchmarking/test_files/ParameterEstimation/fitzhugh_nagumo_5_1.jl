using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "fitzhugh_nagumo_5"
@parameters g a b
@variables t V(t) R(t) y1(t)
D = Differential(t)
states = [V, R]
parameters = [g, a, b]
@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ V,
]
ic = [0.517, 0.432]
p_true = [0.612, 0.215, 0.856]

p_constraints = Dict((g=>(0.0, 1.0)), (a=>(0.0, 1.0)), (b=>(0.0, 1.0)))
ic_constraints = Dict((V=>(0.0, 1.0)), (R=>(0.0, 1.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/fitzhugh_nagumo_5.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
