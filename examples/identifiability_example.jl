import ParameterEstimation

using ModelingToolkit, DifferentialEquations, Plots
using Nemo, HomotopyContinuation
using SIAN: get_order_var
using TaylorSeries

@parameters a01 a21 a12
@variables t x0(t) x1(t) y1(t)
D = Differential(t)

eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]

de = ODESystem(eqs, t, name = :Test)
measured_data = [y1 ~ x0]
identifiability_result = ParameterEstimation.get_identifiability(de;
                                                                 measured_quantities = measured_data)