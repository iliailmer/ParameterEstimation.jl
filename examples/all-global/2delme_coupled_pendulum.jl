using ModelingToolkit, DifferentialEquations

function main()
	@variables xˍt_1 x_1 x_2 x_3 xˍt_2 xˍt_3 alpha
	loss_f =
		(2.0 - x_1)^2 + (3.6 - x_2)^2 + (8.9 - x_3)^2 + 
		(1.07 - xˍt_1)^2 + (1.7 - xˍt_2)^2 + 
		(5.0 - xˍt_3)^2 + (xˍt_1 - alpha * x_1)^2 +
		(xˍt_2 - alpha * x_2)^2 +
		(xˍt_3 - alpha * x_3)^2


	lossvars = get_variables(loss_f)
	f_expr = build_function(loss_f, lossvars)
	n = length(lossvars)
	d = zeros(n)
	g = eval(f_expr)
	g(d)
end

main()

#res = ParameterEstimation.estimate(model, measured_quantities, data_sample, solver=solver)



