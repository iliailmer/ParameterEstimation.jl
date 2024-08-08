using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, ForwardDiff, Plots
using Distributions, Random, StaticArrays
using JLD2, FileIO
solver = Vern9()

@parameters lm d beta a k uu c q b h
@variables t x(t) yy(t) vv(t) w(t) z(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
# TODO
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

ic = [0.388, 0.45, 0.658, 0.148, 0.633]
p_true = [0.312, 0.719, 0.465, 0.555, 0.115, 0.594, 0.59, 0.594, 0.855, 0.645]
time_interval = [-0.5, 0.5]
datasize = 21

sampling_times = range(time_interval[1], time_interval[2], length = datasize)

#prob_true = ODEProblem{false}(model, ic, time_interval, p_true)
#solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times;
#                      abstol = 1e-13, reltol = 1e-13)
#
#data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)
data_sample = load("../data/julia/hiv_2.jld2", "data")

p_rand = rand(Uniform(0.0, 3.0), length(ic) + length(p_true)) # Random Parameters
prob = ODEProblem{false}(model, ic, time_interval, p_rand)
sol = solve(remake(prob, u0 = p_rand[1:length(ic)]), solver,
            p = p_rand[(length(ic) + 1):end],
            saveat = sampling_times;
            abstol = 1e-13, reltol = 1e-13)

function loss(p)
    sol = solve(remake(prob; u0 = SVector{ 5 }(p[1:length(ic)])), Tsit5(), p = p[(length(ic) + 1):end],
                saveat = sampling_times;
                abstol = 1e-13, reltol = 1e-13)
    data_true = [data_sample[v.rhs] for v in measured_quantities]
    data = [(sol[4, :]), (sol[5, :]), (sol[1, :]), (sol[2, :] .+ sol[3, :])]
    if sol.retcode == ReturnCode.Success
        loss = sum(sum((data[i] .- data_true[i]) .^ 2) for i in eachindex(data))
        return loss, sol
    else
        return Inf, sol
    end
end

callback = function (p, l, pred)
    return false
end

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, p_rand, lb = 0.0*ones(5+10), ub = 3.0*ones(5+10))

@time result_ode = Optimization.solve(optprob, BFGS(), callback = callback, maxiters = 200000)

println(result_ode.u)

 all_params = vcat(ic, p_true)
 println("Max. relative abs. error between true and estimated parameters:",
         maximum(abs.((result_ode.u .- all_params) ./ (all_params))))
