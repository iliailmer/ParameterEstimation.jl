using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "crauste_1"
@parameters muN muEE muLE muLL muM muP muPE muPL deltaNE deltaEL deltaLM rhoE rhoP
@variables t n(t) e(t) s(t) m(t) p(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [n, e, s, m, p]
parameters = [muN, muEE, muLE, muLL, muM, muP, muPE, muPL, deltaNE, deltaEL, deltaLM, rhoE, rhoP]
@named model = ODESystem([
                             D(n) ~ -1 * n * muN - n * p * deltaNE,
                             D(e) ~ n * p * deltaNE - e * e * muEE - e * deltaEL + e * p * rhoE,
                             D(s) ~ s * deltaEL - s * deltaLM - s * s * muLL - e * s * muLE,
                             D(m) ~ s * deltaLM - muM * m,
                             D(p) ~ p * p * rhoP - p * muP - e * p * muPE - s * p * muPL,
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ n,
        y2 ~ e,
        y3 ~ s+m,
        y4 ~ p,
]
ic = [0.719, 0.465, 0.555, 0.115, 0.594]
p_true = [0.723, 0.796, 0.883, 0.739, 0.469, 0.724, 0.195, 0.612, 0.215, 0.856, 0.517, 0.432, 0.312]

p_constraints = Dict((muN=>(0.0, 2.0)), (muEE=>(0.0, 2.0)), (muLE=>(0.0, 2.0)), (muLL=>(0.0, 2.0)), (muM=>(0.0, 2.0)), (muP=>(0.0, 2.0)), (muPE=>(0.0, 2.0)), (muPL=>(0.0, 2.0)), (deltaNE=>(0.0, 2.0)), (deltaEL=>(0.0, 2.0)), (deltaLM=>(0.0, 2.0)), (rhoE=>(0.0, 2.0)), (rhoP=>(0.0, 2.0)))
ic_constraints = Dict((n=>(0.0, 2.0)), (e=>(0.0, 2.0)), (s=>(0.0, 2.0)), (m=>(0.0, 2.0)), (p=>(0.0, 2.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/crauste_1.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
