using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations
using Distributions, Random
solver = Tsit5()

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = (0.0, 1.0)
datasize = 100
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

# try HomotopyContinuation for linear algebra

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
                         states, parameters)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = ModelingToolkit.solve(prob_true, solver, p = p_true, saveat = tsteps;
                                      abstol = 1e-12, reltol = 1e-12)

# Initial Parameter Vector
p = abs.(randn!(MersenneTwister(123), zeros(length(p_true) + length(ic))))
_params = Flux.params(p)
prob = ODEProblem(model, [p[1], p[2]], time_interval, p[3:end])

function predict_rd() # Our 1-layer "neural network"
    solve(remake(prob, u0 = [p[1], p[2]]), p = p[3:end], saveat = tsteps)[1, :] # override with new parameters
end

loss_rd() = sum(abs, x - x_true for (x, x_true) in zip(predict_rd(), solution_true[1, :])) # loss function
data = Iterators.repeated((), 3000)
opt = ADAM(0.1)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), AutoTsit5(Rosenbrock23()), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)

# plot(solution_true[1, :], label = "True Solution, r(t)")
# plot!(solution_true[2, :], label = "True Solution, w(t)")

# plot!(predict_rd(), label = "Predicted Solution, r(t)")
# plot!(solve(remake(prob, u0 = [p[1], p[2]]), p = p[3:end], saveat = tsteps)[2, :],
#       label = "Predicted Solution, w(t)")

print("Final Parameter Values: ", p)
print("Mean Abs Error: ", loss_rd())