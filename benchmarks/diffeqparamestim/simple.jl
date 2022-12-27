using ModelingToolkit, Flux, DiffeqParamEstim, DifferentialEquations
using Distributions, Random
solver = Tsit5()

@parameters a b
@variables t x1(t) x2(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2]
parameters = [a, b]

@named model = ODESystem([
                             D(x1) ~ -a * x2,
                             D(x2) ~ 1 / b * (x1),
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ x2,
]

ic = [1.0, 1.0]
p_true = [9.8, 1.3]
time_interval = [0.0, 2.0 * pi * sqrt(1.3 / 9.8)]
datasize = 20
tsteps = range(time_interval[1], time_interval[2], length = datasize)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

data_sample = Dict(v.rhs => solution_true[v.rhs] for v in measured_quantities)

p_rand = rand(Uniform(0.5, 1.5), length(ic) + length(p_true)) # Random Parameters
prob = ODEProblem(model, ic, time_interval,
                  p_rand)
sol = solve(remake(prob, u0 = p_rand[1:length(ic)]), solver,
            p = p_rand[(length(ic) + 1):end],
            saveat = tsteps)

function loss(p)
    sol = solve(remake(prob; u0 = p[1:length(ic)]), Tsit5(), p = p[(length(ic) + 1):end],
                saveat = tsteps)
    data_true = [data_sample[v.rhs] for v in measured_quantities]
    data = [vcat(sol[1, :]) + vcat(sol[3, :]), vcat(sol[2, :])]
    loss = sum(sum((data[i] .- data_true[i]) .^ 2) for i in eachindex(data))
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

result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback, maxiters = 1000)

println(result_ode.u)

all_params = vcat(ic, p_true)
println("Max. relative abs. error between true and estimated parameters:",
        maximum(abs.((result_ode.u .- all_params) ./ (all_params))))