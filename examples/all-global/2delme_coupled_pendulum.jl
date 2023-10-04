using ModelingToolkit, DifferentialEquations
using Optimization
using OptimizationOptimJL
using Zygote

function main()
	@variables xˍt_1 x_1 x_2 x_3 xˍt_2 xˍt_3 alpha
	loss_f =
		(2.0 - x_1)^2 + (3.6 - x_2)^2 + (8.9 - x_3)^2 +
		(1.07 - xˍt_1)^2 + (1.7 - xˍt_2)^2 +
		(5.0 - xˍt_3)^2 + (xˍt_1 - alpha * x_1)^2 +
		(xˍt_2 - alpha * x_2)^2 +
		(xˍt_3 - alpha * x_3)^2


	lossvars = get_variables(loss_f)
	f_expr = build_function(loss_f, lossvars, expression = Val{false})

	#################################  THIS WORKS
	f_expr2(u, p) = f_expr(u)
	u0map = zeros((length(lossvars)))
	g = OptimizationFunction(f_expr2, AutoZygote())
	prob = OptimizationProblem(g, u0map)
	sol = solve(prob, BFGS())
	println("First Version solution:")
	println(sol)
	println(sol.original)
	#########################################3

	#####################################THIS DOESN'T WORK
	@named sys = OptimizationSystem(loss_f, lossvars,[])

	u0dict = Dict()
	for i in lossvars
		u0dict[i] = 0.0
	end

	pnull = Dict()
	prob2=OptimizationProblem(sys, u0dict,pnull)
	sol2 = solve(prob2, NelderMead())

	println("Second Version solution:")
	println(sol2)
	println(sol2.original)
##############################



end

main()

#res = ParameterEstimation.estimate(model, measured_quantities, data_sample, solver=solver)



