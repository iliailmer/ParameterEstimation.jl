using ModelingToolkit, Flux, DiffeqParamEstim, DifferentialEquations
using Distributions, Random
solver = Tsit5()

@parameters lm d beta a k u c q b h
@variables t x(t) y(t) v(t) w(t) z(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [x, y, v, w, z]
parameters = [lm, d, beta, a, k, u, c, q, b, h]

@named model = ODESystem([
                             D(x) ~ lm - d * x - beta * x * v,
                             D(y) ~ beta * x * v - a * y,
                             D(v) ~ k * y - u * v,
                             D(w) ~ c * x * y * w - c * q * y * w - b * w,
                             D(z) ~ c * q * y * w - h * z,
                         ], t, states, parameters)
measured_quantities = [y1 ~ w, y2 ~ z, y3 ~ x, y4 ~ y + v]

ic = [1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [0.0, 10.0]
datasize = 100
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps; abstol = 1e-12,
                      reltol = 1e-12)

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