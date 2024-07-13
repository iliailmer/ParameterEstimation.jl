using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BenchmarkTools
using JLD2, FileIO

solver = Tsit5()

name = "hiv_3"
@parameters lm d beta a k uu c q b h
@variables t x(t) yy(t) vv(t) w(t) z(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [x, yy, vv, w, z]
parameters = [lm, d, beta, a, k, uu, c, q, b, h]
@named model = ODESystem([
                             D(x) ~ lm - d * x - beta * x * vv,
                             D(yy) ~ beta * x * vv - a * yy,
                             D(vv) ~ k * yy - uu * vv,
                             D(w) ~ c * x * yy * w - c * q * yy * w - b * w,
                             D(z) ~ c * q * yy * w - h * z,
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ w,
        y2 ~ z,
        y3 ~ x,
        y4 ~ yy+vv,
]
ic = [0.229, 0.622, 0.303, 0.473, 0.296]
p_true = [0.637, 0.268, 0.203, 0.352, 0.391, 0.556, 0.451, 0.891, 0.182, 0.267]

p_constraints = Dict((lm=>(0.0, 3.0)), (d=>(0.0, 3.0)), (beta=>(0.0, 3.0)), (a=>(0.0, 3.0)), (k=>(0.0, 3.0)), (uu=>(0.0, 3.0)), (c=>(0.0, 3.0)), (q=>(0.0, 3.0)), (b=>(0.0, 3.0)), (h=>(0.0, 3.0)))
ic_constraints = Dict((x=>(0.0, 3.0)), (yy=>(0.0, 3.0)), (vv=>(0.0, 3.0)), (w=>(0.0, 3.0)), (z=>(0.0, 3.0)))
 
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = load("../data/julia/hiv_3.jld2", "data")

#data_sample = ParameterEstimation.sample_data(model, measured_quantities, convert(Array{Float64},time_interval),
#                                              p_true, ic, datasize; solver = solver)

@time res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
            solver = solver, interpolators = Dict("AAA" => ParameterEstimation.aaad), parameter_constraints = p_constraints, ic_constraints = ic_constraints)

all_params = vcat(ic, p_true)
for each in res
  estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
  println("For model ", name, ": Max abs rel. err: ", maximum(abs.((estimates .- all_params) ./ (all_params))))
end
