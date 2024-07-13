using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, ForwardDiff, Plots
using Distributions, Random, StaticArrays
using JLD2, FileIO
solver = Vern9()

@parameters muN muEE muLE muLL muM muP muPE muPL deltaNE deltaEL deltaLM rhoE rhoP
@variables t n(t) e(t) s(t) m(t) p(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
# TODO
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

ic = [0.445, 0.817, 0.394, 0.449, 0.814]
p_true = [0.115, 0.341, 0.628, 0.332, 0.594, 0.443, 0.208, 0.339, 0.556, 0.573, 0.559, 0.623, 0.622]
time_interval = [-0.5, 0.5]
datasize = 21

sampling_times = range(time_interval[1], time_interval[2], length = datasize)

#prob_true = ODEProblem{false}(model, ic, time_interval, p_true)
#solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times;
#                      abstol = 1e-13, reltol = 1e-13)
#
#data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)
data_sample = load("../data/julia/crauste_7.jld2", "data")

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
    data = [(sol[1, :]), (sol[2, :]), (sol[3, :] .+ sol[4, :]), (sol[5, :])]
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
optprob = Optimization.OptimizationProblem(optf, p_rand, lb = 0.0*ones(5+13), ub = 3.0*ones(5+13))

@time result_ode = Optimization.solve(optprob, BFGS(), callback = callback, maxiters = 200000)

println(result_ode.u)

 all_params = vcat(ic, p_true)
 println("Max. relative abs. error between true and estimated parameters:",
         maximum(abs.((result_ode.u .- all_params) ./ (all_params))))
