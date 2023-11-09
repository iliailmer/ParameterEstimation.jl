module ParameterEstimation

using PrecompileTools
@recompile_invalidations begin

	import DifferentialEquations: Tsit5
	import TaylorSeries: Taylor1
	import OrderedCollections: OrderedDict
	using BaryRational: BaryRational
	using ForwardDiff: ForwardDiff
	using LinearAlgebra: LinearAlgebra
	using TaylorDiff: TaylorDiff
	using Suppressor

	using ProgressMeter, Logging, Printf
	using ModelingToolkit, LinearSolve, LinearAlgebra
	using SIAN, HomotopyContinuation, Groebner, Oscar
	using .ReturnCode
	import StructuralIdentifiability: eval_at_nemo, ODE
	using BaryRational
	using ForwardDiff
	using ArbNumerics

end
Float = Union{Float64, Float32, Float16}
include("includes.jl")

export check_identifiability, estimate, filter_solutions
@recompile_invalidations begin
	@compile_workload begin

		@parameters a b c d
		@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
		D = Differential(t)
		states = [x1, x2, x3, x4]
		parameters = [a, b, c, d]

		@named model = ODESystem([
				D(x1) ~ a + x2,
				D(x2) ~ b + x3,
				D(x3) ~ c + x4,
				D(x4) ~ d,
			], t, states, parameters)
		measured_quantities = [
			y1 ~ x1,
			y2 ~ x2,
			y3 ~ x3,
			y4 ~ x4,
		]

		ic = [0.0, 0.0, 0.0, 0.0]
		p_true = [2.0, 3.0, 4.0, 5.0]
		time_interval = [-4.0, 4.0]
		datasize = 9
		data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
			p_true, ic, datasize;)
		res = ParameterEstimation.estimate(model, measured_quantities, data_sample)


	end
end
end

