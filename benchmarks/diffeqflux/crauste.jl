using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = Tsit5()

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

ic = [1.0, -1.0, 1.0, -1.0, 1.0]
time_interval = [0.0, 1.0]
datasize = 10
tsteps = range(time_interval[1], time_interval[2], length = datasize)

p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5, 1.0, 1.0, 1.0, 1.0, 0.9, 1.2] # True Parameters
prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2], p[3], p[4], p[5]], time_interval, p[6:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2], p[3], p[4], p[5]]), solver, p = p[6:end],
                saveat = tsteps) # override with new parameters
    return [sol[1, :] sol[2, :] sol[3, :] + sol[4, :] sol[5, :]]
end

function loss_rd()
    y_pred = predict_rd()
    y_true = [solution_true[1, :] solution_true[2, :] solution_true[3, :] +
                                                      solution_true[4, :] solution_true[5,
                                                                                        :]]
    return sum(abs, y_true .- y_pred)
end # loss function
data = Iterators.repeated((), 3000)
opt = ADAM(0.01)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), AutoTsit5(Rosenbrock23()), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)
println("parameters:", p)
# plot(solution_true[1, :], label = "True Solution")
# plot!(predict_rd()[:, 1], label = "Predicted Solution")

# plot!(solution_true[2, :], label = "True Solution")
# plot!(predict_rd()[:, 2], label = "Predicted Solution")

# plot!(solution_true[3, :] + solution_true[4, :], label = "True Solution")
# plot!(predict_rd()[:, 3], label = "Predicted Solution")

# plot!(solution_true[5, :], label = "True Solution")
# plot!(predict_rd()[:, 4], label = "Predicted Solution")