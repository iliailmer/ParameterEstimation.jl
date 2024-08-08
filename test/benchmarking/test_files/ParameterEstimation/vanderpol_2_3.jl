using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "vanderpol_2"
@parameters a b
@variables t x1(t) x2(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2]
parameters = [a, b]
@named model = ODESystem([
                             D(x1) ~ a * x2,
                             D(x2) ~ -(x1) - b * (x1^2 - 1) * (x2),
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ x1,
        y2 ~ x2,
]
ic = [0.733, 0.523]
p_true = [0.871, 0.407]

p_constraints = Dict((a=>(0.0, 3.0)), (b=>(0.0, 3.0)))
ic_constraints = Dict((x1=>(0.0, 3.0)), (x2=>(0.0, 3.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/vanderpol_2.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
