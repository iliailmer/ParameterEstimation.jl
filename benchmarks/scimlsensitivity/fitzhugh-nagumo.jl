using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, ForwardDiff, Plots
using Distributions, Random, StaticArrays
solver = Vern9()

@parameters g a b
@variables t V(t) R(t) y1(t) y2(t)
D = Differential(t)
states = [V, R]
parameters = [g, a, b]

ic = SA[1.0, -1.0]
time_interval = [0.0, 1]
datasize = 50
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [2, 2 / 10, 2 / 10] # True Parameters
measured_quantities = [y1 ~ V]

@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ], t, states, parameters)

prob_true = ODEProblem{false}(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times;
                      abstol = 1e-10, reltol = 1e-10)

data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)

p_rand = rand(Uniform(0.5, 1.5), length(ic) + length(p_true)) # Random Parameters
prob = ODEProblem{false}(model, ic, time_interval,
                  p_rand)
sol = solve(remake(prob, u0 = p_rand[1:length(ic)]), solver,
            p = p_rand[(length(ic) + 1):end],
            saveat = sampling_times; abstol = 1e-10, reltol = 1e-10)

function loss(p)
    sol = solve(remake(prob; u0 = SVector{2}(p[1:length(ic)])), Tsit5(), p = p[(length(ic) + 1):end],
                saveat = sampling_times)
    data_true = [data_sample[v.rhs] for v in measured_quantities]
    data = [sol[1, :]]
    if sol.retcode == ReturnCode.Success
        loss = sum(sum((data[i] .- data_true[i]) .^ 2) for i in eachindex(data))
        return loss, sol
    else
        return Inf, sol
    end
end

callback = function (p, l, pred)
    # display(l)
    #     plt = plot(pred, ylim = (0, 6))
    #     display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, p_rand)

result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback, maxiters = 1000)

println(result_ode.u)

all_params = vcat(ic, p_true)
println("Max. relative abs. error between true and estimated parameters:",
        maximum(abs.((result_ode.u .- all_params) ./ (all_params))))

num_unknowns = length(ic) + length(p_true)
all_params = vcat(ic, p_true)
using OrderedCollections
size_err_map = OrderedDict{Int, Float64}()
for datasize in 3:21
    sampling_times = range(time_interval[1], time_interval[2], length = datasize)
    prob_true = ODEProblem(model, ic, time_interval, p_true)
    solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times)
    data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)
    p_rand = rand(Uniform(0.5, 1.5), num_unknowns) # Random Parameters
    prob = ODEProblem(model, ic, time_interval, p_rand)
    optprob = Optimization.OptimizationProblem(optf, p_rand)
    result_ode = p_rand
    try
        result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback,
                                        maxiters = 1000)
    catch
        result_ode = p_rand # in case of failure use random parameters as answer
    end
    size_err_map[datasize] = maximum(100 *
                                     abs.((result_ode .- all_params) ./ (all_params)))
end
open("fhn_t_$(time_interval[1])_$(time_interval[2]).txt", "w") do f
    for (k, v) in size_err_map
        println(f, "$k $v")
    end
end
