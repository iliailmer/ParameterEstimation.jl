using ModelingToolkit
using DifferentialEquations, Plots

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

u0 = [100.0, 100.0]
tspan = (0.0, 10.0)
datasize = 100
tsteps = range(tspan[1], tspan[2], length=datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters

@named lv = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ -k3 * w + k2 * r * w])

prob_true = ODEProblem(lv, u0, tspan, p_true)
solution_true = solve(prob_true, Tsit5(), p=p_true, saveat=tsteps)
