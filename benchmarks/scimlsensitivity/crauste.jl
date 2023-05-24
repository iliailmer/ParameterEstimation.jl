using ModelingToolkit, DifferentialEquations, Optimization, OptimizationPolyalgorithms,
      OptimizationOptimJL, SciMLSensitivity, ForwardDiff, Plots
using Distributions, Random, StaticArrays
solver = Vern9()

@parameters mu_N mu_EE mu_LE mu_LL mu_M mu_P mu_PE mu_PL delta_NE delta_EL delta_LM rho_E rho_P
@variables t N(t) E(t) S(t) M(t) P(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [N, E, S, M, P]
parameters = [
    mu_N,
    mu_EE,
    mu_LE,
    mu_LL,
    mu_M,
    mu_P,
    mu_PE,
    mu_PL,
    delta_NE,
    delta_EL,
    delta_LM,
    rho_E,
    rho_P,
]
@named model = ODESystem([
                             D(N) ~ -N * mu_N - N * P * delta_NE,
                             D(E) ~ N * P * delta_NE - E^2 * mu_EE -
                                    E * delta_EL + E * P * rho_E,
                             D(S) ~ S * delta_EL - S * delta_LM - S^2 * mu_LL -
                                    E * S * mu_LE,
                             D(M) ~ S * delta_LM - mu_M * M,
                             D(P) ~ P^2 * rho_P - P * mu_P - E * P * mu_PE -
                                    S * P * mu_PL,
                         ], t, states, parameters)
measured_quantities = [y1 ~ N, y2 ~ E, y3 ~ S + M, y4 ~ P]

ic = SA[1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [0.0, 1.0]
datasize = 20
sampling_times = range(time_interval[1], time_interval[2], length = datasize)

p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5, 1.0, 1.0, 1.0, 1.0, 0.9, 1.2] # True Parameters
prob_true = ODEProblem{false}(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times;
                      abstol = 1e-10, reltol = 1e-10)
data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)

p_rand = rand(Uniform(0.5, 1.5), length(ic) + length(p_true)) # Random Parameters
prob = ODEProblem{false}(model, ic, time_interval,
                  p_rand)
# sol = solve(remake(prob, u0 = p[1:length(ic)]), solver, p = p_rand[(length(ic) + 1):end],
# saveat = sampling_times)

function loss(p)
    sol = solve(remake(prob; u0 = SVector{5}(p[1:length(ic)])), Tsit5(), p = p[(length(ic) + 1):end],
                saveat = sampling_times;
                abstol = 1e-10, reltol = 1e-10)
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

@time result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback, maxiters = 1000)

println(result_ode.u)

all_params = vcat(ic, p_true)
println("Max. relative abs. error between true and estimated parameters:",
        maximum(abs.((result_ode.u .- all_params) ./ (all_params))))
# num_unknowns = length(ic) + length(p_true)
# all_params = vcat(ic, p_true)
# using OrderedCollections
# size_err_map = OrderedDict{Int, Float64}()
# for datasize in 3:21
#     prob_true = ODEProblem(model, ic, time_interval, p_true)
#     solution_true = solve(prob_true, solver, p = p_true, saveat = sampling_times)
#     data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)
#     p_rand = rand(Uniform(0.5, 1.5), length(ic) + length(p_true)) # Random Parameters
#     prob = ODEProblem(model, ic, time_interval,
#                       p_rand)
#     optprob = Optimization.OptimizationProblem(optf, p_rand)

#     result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback,
#                                     maxiters = 1000)

#     size_err_map[datasize] = maximum(100 *
#                                      abs.((result_ode .- all_params) ./ (all_params)))
# end

# using Plots
# scatter(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
#         title = "Crauste Model, $(num_unknowns) unknowns", legend = false)
# plot!(size_err_map, xlabel = "Number of data points", ylabel = "Max. rel. err. [%]",
#       title = "Crauste Model, $(num_unknowns) unknowns", legend = false)
# # save size_err_map to file
# open("Crauste_t_$(time_interval[1])_$(time_interval[2]).txt", "w") do f
#     for (k, v) in size_err_map
#         println(f, "$k $v")
#     end
# end
# # png("Crauste_t_$(time_interval[1])_$(time_interval[2]).png")
