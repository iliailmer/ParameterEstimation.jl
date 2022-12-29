using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, Zygote, Plots
using Distributions, Random
solver = Tsit5()

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = (0.0, 1.0)
datasize = 8
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_quantities = [r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
                         states, parameters)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, solver, p = p_true, saveat = tsteps;
                                      abstol = 1e-12, reltol = 1e-12)
data_sample = [solution_true[v] for v in measured_quantities]

p_rand = rand(Uniform(0.02, 1.5), length(p_true)) # Random Parameters
prob = ODEProblem(model, ic, time_interval,
                  p_rand)
sol = solve(prob, solver, p = p_rand, saveat = tsteps)

function loss(p)
    sol = solve(prob, Tsit5(), p = p, saveat = tsteps)
    data = [sol[v] for v in measured_quantities]
    loss = sum(sum((data[i] .- data_sample[i]) .^ 2 for i in eachindex(data)))
    return loss, sol
end

callback = function (p, l, pred)
    display(l)
    #     plt = plot(pred, ylim = (0, 6))
    #     display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, p_rand)

result_ode = Optimization.solve(optprob, PolyOpt(),
                                callback = callback,
                                maxiters = 100)