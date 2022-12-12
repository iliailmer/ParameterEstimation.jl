using Pkg
Pkg.activate(; temp = true)
Pkg.add(["ModelingToolkit", "Flux", "DiffEqFlux", "DifferentialEquations", "Plots"])
Pkg.add(["Distributions", "Random"])

using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = AutoTsit5(Rosenbrock23())

@parameters mu bi bw al g dz k
@variables t s(t) i(t) w(t) r(t) y1(t) y2(t)
D = Differential(t)

states = [s, i, w, r]
parameters = [mu, bi, bw, al, g, dz, k]

@named model = ODESystem([
                             D(s) ~ mu - bi * s * i - bw * s * w - mu * s + al * r,
                             D(i) ~ bw * s * w + bi * s * i - g * i - mu * i,
                             D(w) ~ dz * (i - w),
                             D(r) ~ g * i - mu * r - al * r,
                         ], t, states, parameters)
measured_quantities = [y1 ~ i, y2 ~ i + r + s]

ic = [1.0, -1.0, 1.0, -1.0]
time_interval = [0.0, 1.0]
datasize = 10
tsteps = range(time_interval[1], time_interval[2], length = datasize)

p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5] # True Parameters
prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2], p[3], p[4]], time_interval, p[5:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2], p[3], p[4]]), solver, p = p[5:end],
                saveat = tsteps) # override with new parameters
    return [sol[2, :] sol[2, :] + sol[4, :] + sol[1, :]]
end

function loss_rd()
    pred = predict_rd()
    y_true = [solution_true[2, :] solution_true[2, :] + solution_true[4, :] +
                                  solution_true[1, :]]
    return sum(abs2, y_true .- pred)
end # loss function
data = Iterators.repeated((), 2500)
opt = ADAM(0.1)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), Tsit5(), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)

plot(solution_true[1, :], label = "True Solution")
plot!(predict_rd()[:, 1], label = "Predicted Solution")

plot(solution_true[2, :] + solution_true[4, :] + solution_true[1, :],
     label = "True Solution")
plot!(predict_rd()[:, 2], label = "Predicted Solution")
