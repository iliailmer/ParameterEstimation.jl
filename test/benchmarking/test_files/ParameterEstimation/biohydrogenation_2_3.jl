using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "biohydrogenation_2"
@parameters k5 k6 k7 k8 k9 k10
@variables t x4(t) x5(t) x6(t) x7(t) y1(t) y2(t)
D = Differential(t)
states = [x4, x5, x6, x7]
parameters = [k5, k6, k7, k8, k9, k10]
@named model = ODESystem([
                             D(x4) ~ - k5 * x4 / (k6 + x4),
                             D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5/(k8 + x5 + x6),
                             D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
                             D(x7) ~ k9 * x6 * (k10 - x6) / k10,
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ x4,
        y2 ~ x5,
]
ic = [0.215, 0.856, 0.517, 0.432]
p_true = [0.883, 0.739, 0.469, 0.724, 0.195, 0.612]

p_constraints = Dict((k5=>(0.0, 3.0)), (k6=>(0.0, 3.0)), (k7=>(0.0, 3.0)), (k8=>(0.0, 3.0)), (k9=>(0.0, 3.0)), (k10=>(0.0, 3.0)))
ic_constraints = Dict((x4=>(0.0, 3.0)), (x5=>(0.0, 3.0)), (x6=>(0.0, 3.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/biohydrogenation_2.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
