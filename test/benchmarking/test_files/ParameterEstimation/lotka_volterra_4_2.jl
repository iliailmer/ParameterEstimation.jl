using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "lotka_volterra_4"
@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)
states = [r, w]
parameters = [k1, k2, k3]
@named model = ODESystem([
                             D(r) ~ k1*r - k2*r*w,
                             D(w) ~ k2*r*w - k3*w,
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ r,
]
ic = [0.724, 0.195]
p_true = [0.883, 0.739, 0.469]

p_constraints = Dict((k1=>(0.0, 2.0)), (k2=>(0.0, 2.0)), (k3=>(0.0, 2.0)))
ic_constraints = Dict((r=>(0.0, 2.0)), (w=>(0.0, 2.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/lotka_volterra_4.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
