using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, ForwardDiff, Plots
using Distributions, Random, StaticArrays
using JLD2, FileIO
solver = Vern9()

@parameters a b nu
@variables t S(t) E(t) In(t) NN(t) y1(t) y2(t)
D = Differential(t)
# TODO
states = [S, E, In, NN]
parameters = [a, b, nu]
@named model = ODESystem([
                             D(S) ~ -b * S * In / NN,
                             D(E) ~ b * S * In / NN - nu * E,
                             D(In) ~ nu * E - a * In,
                             D(NN) ~ 0,
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ In,
        y2 ~ NN,
]

ic = [0.766, 0.723, 0.796, 0.883]
p_true = [0.157, 0.17, 0.116]
time_interval = [-0.5, 0.5]
datasize = 21

sampling_times = range(time_interval[1], time_interval[2], length = datasize)

#prob_true = ODEProblem{false}(model, ic, time_interval, p_true)
#solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times;
#                      abstol = 1e-13, reltol = 1e-13)
#
#data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)
data_sample = load("../data/julia/seir_2.jld2", "data")

p_rand = rand(Uniform(0.0, 2.0), length(ic) + length(p_true)) # Random Parameters
prob = ODEProblem{false}(model, ic, time_interval, p_rand)
sol = solve(remake(prob, u0 = p_rand[1:length(ic)]), solver,
            p = p_rand[(length(ic) + 1):end],
            saveat = sampling_times;
            abstol = 1e-13, reltol = 1e-13)

function loss(p)
    sol = solve(remake(prob; u0 = SVector{ 4 }(p[1:length(ic)])), Tsit5(), p = p[(length(ic) + 1):end],
                saveat = sampling_times;
                abstol = 1e-13, reltol = 1e-13)
    data_true = [data_sample[v.rhs] for v in measured_quantities]
    data = [(sol[3, :]), (sol[4, :])]
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
optprob = Optimization.OptimizationProblem(optf, p_rand, lb = 0.0*ones(4+3), ub = 2.0*ones(4+3))

@time result_ode = Optimization.solve(optprob, BFGS(), callback = callback, maxiters = 200000)

println(result_ode.u)

 all_params = vcat(ic, p_true)
 println("Max. relative abs. error between true and estimated parameters:",
         maximum(abs.((result_ode.u .- all_params) ./ (all_params))))
