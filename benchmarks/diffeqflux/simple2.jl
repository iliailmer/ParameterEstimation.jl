using ModelingToolkit, DifferentialEquations, Plots

using ParameterEstimation

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
datasize = 10
tsteps = range(time_interval[1], time_interval[2], length = datasize)

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
data = Iterators.repeated((), 1000)
opt = ADAM(0.01)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), Tsit5(), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)

plot(solution_true[1, :], label = "True: x1")
plot!(solution_true[2, :], label = "True: x2")
plot!(predict_rd()[:, 1], label = "Predicted: x1")
plot!(predict_rd()[:, 2], label = "Predicted: x2")