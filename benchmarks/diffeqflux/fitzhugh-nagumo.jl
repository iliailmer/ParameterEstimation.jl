using Pkg
Pkg.activate(; temp = true)
Pkg.add(["ModelingToolkit", "Flux", "DiffEqFlux", "DifferentialEquations", "Plots"])
Pkg.add(["Distributions", "Random"])

using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = AutoTsit5(Rosenbrock23())

function fhn(du, u, p, t)
    du[1] = p[1] * (u[1] - u[1]^3 / 3 + u[2])
    du[2] = -1 / p[1] * (u[1] - p[2] + p[3] * u[2])
end

@parameters g a b
@variables t V(t) R(t) y1(t) y2(t)
D = Differential(t)
states = [V, R]
parameters = [g, a, b]

ic = [1.0, -1.0]
time_interval = (0.0, 5)
datasize = 50
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [2, 2 / 10, 2 / 10] # True Parameters
measured_quantities = [y1 ~ V]

@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ], t, states, parameters)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)
prob = ODEProblem(model, [p[1], p[2]], time_interval, p[3:end])

function predict_rd() # Our 1-layer "neural network"
    solve(remake(prob, u0 = [p[1], p[2]]), solver, p = p[3:end], saveat = tsteps)[1, :] # override with new parameters
end

loss_rd() = sum(abs2, x - x_true for (x, x_true) in zip(predict_rd(), solution_true[1, :])) # loss function
data = Iterators.repeated((), 2500)
opt = ADAM(0.1)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), AutoTsit5(Rosenbrock23()), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)

plot(solution_true[1, :], label = "True Solution")
plot!(predict_rd(), label = "Predicted Solution")
