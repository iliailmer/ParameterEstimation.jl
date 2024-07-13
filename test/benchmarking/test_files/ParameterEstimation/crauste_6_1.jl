using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "crauste_6"
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
ic = [0.501, 0.865, 0.615, 0.439, 0.585]
p_true = [0.278, 0.862, 0.458, 0.777, 0.66, 0.338, 0.751, 0.417, 0.805, 0.565, 0.805, 0.654, 0.68]

p_constraints = Dict((muN=>(0.0, 1.0)), (muEE=>(0.0, 1.0)), (muLE=>(0.0, 1.0)), (muLL=>(0.0, 1.0)), (muM=>(0.0, 1.0)), (muP=>(0.0, 1.0)), (muPE=>(0.0, 1.0)), (muPL=>(0.0, 1.0)), (deltaNE=>(0.0, 1.0)), (deltaEL=>(0.0, 1.0)), (deltaLM=>(0.0, 1.0)), (rhoE=>(0.0, 1.0)), (rhoP=>(0.0, 1.0)))
ic_constraints = Dict((n=>(0.0, 1.0)), (e=>(0.0, 1.0)), (s=>(0.0, 1.0)), (m=>(0.0, 1.0)), (p=>(0.0, 1.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/crauste_6.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
