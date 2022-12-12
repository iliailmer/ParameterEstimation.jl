using Pkg
Pkg.activate(; temp = true)
Pkg.add(["ModelingToolkit", "Flux", "DiffEqFlux", "DifferentialEquations", "Plots"])
Pkg.add(["Distributions", "Random"])

using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = AutoTsit5(Rosenbrock23())

@parameters a b c d
@variables t x1(t) x2(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2]
parameters = [a, b, c, d]

@named model = ODESystem([
                             D(x1) ~ a * x1 + b * x2,
                             D(x2) ~ c * x1 + d * x2,
                         ], t, states, parameters)
measured_quantities = [
    y1 ~ x1,
    y2 ~ x2,
]

ic = [1.0, 1.0]
time_interval = (0.0, 10.0)
datasize = 10
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [1.5, -2, 1, -1.0]

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2]], time_interval, p[3:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2]]), solver,
                p = p[3:end], saveat = tsteps) # override with new parameters
    return [sol[1, :] sol[2, :]]
end

function loss_rd()
    y_pred = predict_rd()
    y_true = [solution_true[1, :] solution_true[2, :]]
    return sum(abs2, y_true .- y_pred)
end # loss function
data = Iterators.repeated((), 2500)
opt = ADAM(0.0001)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), Tsit5(), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)
