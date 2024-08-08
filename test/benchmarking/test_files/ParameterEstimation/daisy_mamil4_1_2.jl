using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "daisy_mamil4_1"
@parameters k01 k12 k13 k14 k21 k31 k41
@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t)
D = Differential(t)
states = [x1, x2, x3, x4]
parameters = [k01, k12, k13, k14, k21, k31, k41]
@named model = ODESystem([
                             D(x1) ~ -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 - k41 * x1,
                             D(x2) ~ -k12 * x2 + k21 * x1,
                             D(x3) ~ -k13 * x3 + k31 * x1,
                             D(x4) ~ -k14 * x4 + k41 * x1,
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ x1,
        y2 ~ x2,
        y3 ~ x3 + x4,
]
ic = [0.723, 0.796, 0.883, 0.739]
p_true = [0.523, 0.554, 0.84, 0.157, 0.17, 0.116, 0.766]

p_constraints = Dict((k01=>(0.0, 2.0)), (k12=>(0.0, 2.0)), (k13=>(0.0, 2.0)), (k14=>(0.0, 2.0)), (k21=>(0.0, 2.0)), (k31=>(0.0, 2.0)), (k41=>(0.0, 2.0)))
ic_constraints = Dict((x1=>(0.0, 2.0)), (x2=>(0.0, 2.0)), (x3=>(0.0, 2.0)), (x4=>(0.0, 2.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/daisy_mamil4_1.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
