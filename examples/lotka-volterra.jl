import ParameterEstimation

using ModelingToolkit, DifferentialEquations, Plots
using Nemo, HomotopyContinuation
using SIAN: get_order_var
using TaylorSeries

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

u0 = [100.0, 100.0]
tspan = (0.0, 1.0)
datasize = 20
tsteps = range(tspan[1], tspan[2], length=datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_data = [y1 ~ r]

@named lv = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ -k3 * w + k2 * r * w])

prob_true = ODEProblem(lv, u0, tspan, p_true)
solution_true = ModelingToolkit.solve(prob_true, Tsit5(), p=p_true, saveat=tsteps)
r_sample = solution_true[1, :]
g = ParameterEstimation.interpolate(tsteps, r_sample, Int(datasize / 2))
err = maximum(abs.(r_sample - g.(tsteps)))
print("Error in interpolation: ", err,)
diff_order = length(p_true) + length(u0)

τ = Taylor1(diff_order + 1)
g_expanded = g(τ)

id_result = ParameterEstimation.get_identifiability(lv; measured_quantities=measured_data)

# convert indeterminates to symbols in HomotopyContinuation

