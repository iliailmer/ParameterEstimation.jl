using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, ForwardDiff, Plots
using Distributions, Random, StaticArrays
solver = Vern9()

@parameters a12 a13 a21 a31 a01
@variables t x1(t) x2(t) x3(t) y1(t) y2(t)
D = Differential(t)

ic = SA[1.0, 2.0, 1.0]
time_interval = [0.0, 1.0]
datasize = 10
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.2, 0.3, 0.5, 0.6, -0.2] # True Parameters

states = [x1, x2, x3]
parameters = [a12, a13, a21, a31, a01]
@named model = ODESystem([D(x1) ~ -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3,
                             D(x2) ~ a21 * x1 - a12 * x2,
                             D(x3) ~ a31 * x1 - a13 * x3],
                         t, states, parameters)
measured_quantities = [y1 ~ x1, y2 ~ x2]
prob_true = ODEProblem{false}(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times;
                      abstol = 1e-10, reltol = 1e-10)
data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)

p_rand = rand(Uniform(0.5, 1.5), length(ic) + length(p_true)) # Random Parameters
prob = ODEProblem{false}(model, ic, time_interval, p_rand)
sol = solve(remake(prob, u0 = p_rand[1:length(ic)]), solver,
            p = p_rand[(length(ic) + 1):end],
            saveat = sampling_times; abstol = 1e-10, reltol = 1e-10)

function loss(p)
    sol = solve(remake(prob; u0 = SVector{3}(p[1:length(ic)])), Tsit5(), p = p[(length(ic) + 1):end],
                saveat = sampling_times; abstol = 1e-10, reltol = 1e-10)
    data_true = [data_sample[v.rhs] for v in measured_quantities]
    data = [vcat(sol[1, :]), vcat(sol[2, :])]
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
    prob_true = ODEProblem(model, ic, time_interval, p_true)
    solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times)
    data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)
    p_rand = rand(Uniform(0.5, 1.5), length(ic) + length(p_true)) # Random Parameters
    prob = ODEProblem(model, ic, time_interval,
                      p_rand)
    optprob = Optimization.OptimizationProblem(optf, p_rand)

    result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback,
                                    maxiters = 1000)

    size_err_map[datasize] = maximum(100 *
                                     abs.((result_ode .- all_params) ./ (all_params)))
end
open("daisy_mamil3_t_$(time_interval[1])_$(time_interval[2]).txt", "w") do f
    for (k, v) in size_err_map
        println(f, "$k $v")
    end
end
