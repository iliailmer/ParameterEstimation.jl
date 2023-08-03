using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BaryRational
using ForwardDiff
import HomotopyContinuation as HC
using Optimization
using OptimizationOptimJL


######### TODO(orebas)  REFACTOR into bary_derivs or something similiar
#to use the below, you can just pass vectors of xvalues and yvalues like so:
#F = aaad(xdata, ydata)
#and then F will be a callable, i.e. F(0.5) should work.  Let's please restrict to real xdata, and if so it should be sorted.  F is only defined in the range of the xvalues.

#To construct a derivative, you can do
#derivf(z) = ForwardDiff.derivative(F, z)
#I hope and suspect that other AD frameworks should work as well.

function baryEval(z, f::Vector{T}, x::Vector{T}, w::Vector{T}, tol = 1e-13) where {T}
	@assert(length(f) == length(x) == length(w))
	num = zero(T)
	den = zero(T)
	breakflag = false
	breakindex = -1
	for j in eachindex(f)
		t = w[j] / (z - x[j])
		num += t * f[j]
		den += t
		if (abs(z - x[j]) < tol)
			breakflag = true
			breakindex = j
		end
	end
	fz = num / den
	if (breakflag)
		num = zero(T)
		den = zero(T)
		for j in eachindex(f)
			if (j != breakindex)
				t = w[j] / (z - x[j])
				num += t * f[j]
				den += t
			end
		end
		m = z - x[breakindex]
		fz = (w[breakindex] * f[breakindex] + m * num) / (w[breakindex] + m * den)
	end
	return fz
end

struct AAADapprox{T}
	internalAAA::T
end

(y::AAADapprox)(z) = baryEval(z, y.internalAAA.f, y.internalAAA.x, y.internalAAA.w)

function aaad(xs::AbstractArray{T}, ys::AbstractArray{T}) where {T}
	@assert length(xs) == length(ys)
	internalApprox = aaa(xs, ys)
	return AAADapprox(internalApprox)
end

####################





function SimpleParameterEstimation(model::ODESystem, measured_quantities, data_sample, solver)
	println("Starting")
	#build the equation array.  
	#eqns[i,*] relates to the time index i, i.e. 1 is the first time and sample_count is the last time
	#eqns[i,1] is the equations coming from the model, enforced with respect to time index i
	#eqns[i,2] is the sample data we are given with respect to time index i:  all the measured quantities.
	#eqns[i,3] is the additional equations we get from taking first derivatives of the measured quantities, and assuming they are
	#equal to derivatives estimated numerically via the given interpolator.  
	t = ModelingToolkit.get_iv(model)
	model_eq = equations(model)
	model_states = states(model)
	model_ps = parameters(model)


	t_vector = pop!(data_sample, "t") #TODO(orebas) make it use the independent variable name
	sample_count = length(t_vector)
	D = Differential(t)

	eqns = []

	interpolants = Dict()
	for j in measured_quantities  #this loop constructs the interpolators of measured data
		#TODO this seems like a brittle design.  
		#We are relying on an exact match between the measured quantities and the data_sample
		#data_vector = data_sample(j.rhs)
		r = j.rhs
		y_vector = data_sample[r]
		interpolants[r] = aaad(t_vector, y_vector)
	end


	for i in 1:sample_count
		equations_time_slice_full = []

		d = Dict()

		for j in eachindex(model_states)
			varname = (model_states[j].metadata.value[2])
			vname2 = Symbol("$(varname)_$(i)")
			vname3 = @variables $vname2
			d[model_states[j]] = (vname3[1])

			dvar = ModelingToolkit.diff2term(expand_derivatives(D(model_states[j])))
			dvarname = dvar.metadata.value[2]
			dvname2 = Symbol("$(dvarname)_$(i)")
			dvname3 = @variables $dvname2
			d[dvar] = dvname3[1]

		end
		equations_time_slice_from_ODE_only = []# assemple the ODE equations
		equations_time_slice_from_measured_quantities_0th_deriv = []
		equations_time_slice_from_measured_quantities_1st_deriv = []

		for j in eachindex(model_eq)
			lhs1 = expand_derivatives(model_eq[j].lhs)
			rhs1 = expand_derivatives(model_eq[j].rhs)
			lhs2 = ModelingToolkit.diff2term(lhs1)
			rhs2 = ModelingToolkit.diff2term(rhs1)
			lhs3 = substitute(lhs2, d)
			rhs3 = substitute(rhs2, d)
			#push!(equations_time_slice_from_ODE_only, lhs3 ~ rhs3)
			push!(equations_time_slice_from_ODE_only, lhs3 - rhs3)

		end


		for j in measured_quantities
			r = j.rhs
			r1 = ModelingToolkit.diff2term(expand_derivatives(r))
			r2 = substitute(r1, d)
			yval = interpolants[r](t_vector[i])
			#eq = yval ~ r2
			eq = yval - r2

			push!(equations_time_slice_from_measured_quantities_0th_deriv, eq)
		end

		for j in measured_quantities

			r = j.rhs
			dr = D(r)
			dr1 = ModelingToolkit.diff2term(expand_derivatives(dr))
			dr2 = substitute(dr1, d)
			yval = ForwardDiff.derivative(interpolants[r], t_vector[i])
			#eq = yval ~ dr2
			eq = yval - dr2

			push!(equations_time_slice_from_measured_quantities_1st_deriv, eq)
		end
		#println("ODE ONLY")
		#println(equations_time_slice_from_ODE_only)
		#println("0th deriv of measured quantities:")

		#println(equations_time_slice_from_measured_quantities_0th_deriv)
		#println("1st deriv of measured quantites")

		#println(equations_time_slice_from_measured_quantities_1st_deriv)

		push!(equations_time_slice_full, equations_time_slice_from_ODE_only)
		push!(equations_time_slice_full, equations_time_slice_from_measured_quantities_0th_deriv)
		push!(equations_time_slice_full, equations_time_slice_from_measured_quantities_1st_deriv)
		#push!(eqns, model_converted)
		#println(measured_quantities)
		#println(data_sample)
		push!(eqns, equations_time_slice_full)

	end
	for i in eachindex(eqns)
		for j in eachindex(eqns[i])
			println()
			println("$i, $j")
			println(eqns[i][j])
		end
	end
	loss = typeof(eqns[1][1][1])(0)

	for i in eachindex(eqns)
		for j in eachindex(eqns[i])
			for k in eachindex(eqns[i][j])
				loss += (eqns[i][j][k])^2
			end
		end
	end

	println(loss)
	println(typeof(loss))

	lossvars = get_variables(loss)

	println(lossvars)

	f_expr = build_function(loss, lossvars, expression = Val{false})
	f_expr2(u, p) = f_expr(u)


	n = length(lossvars)
u0map = zeros(n)
	

	println(f_expr)
	println(typeof(f_expr))

	println(typeof(loss))
	@syms loss
	println(f_expr(u0map))

	@named syst = OptimizationSystem(f_expr2, lossvars, [])
	u0map = Dict()
	for i in lossvars
		u0map[i] = 0
	end
	u0map = zeros(length(lossvars))
	loss2(u, p) = loss(u)
	println(f_expr(u0map))
	#f = OptimizationFunction(loss2)
	prob = OptimizationProblem(f_expr2, u0map)
	sol = solve(prob, SimulatedAnnealing())
	print(sol)
end

function main()
	solver = Tsit5()

	@parameters alpha
	@variables t theta_1(t) y1(t)
	D = Differential(t)

	#coupled pendulum equations
	states = [theta_1]
	parameters = [alpha]
	@named model = ODESystem([
			D(theta_1) ~ alpha * theta_1], t, states, parameters)


	#initial conditions
	ic = [2.0]
	p_true = [0.5]
	time_interval = [0.0, 3]
	datasize = 4

	v = randn(datasize)
	v = sort((v .- minimum(v)) / (maximum(v) - minimum(v))) * time_interval[2]

	measured_quantities = [y1 ~ theta_1]
	data_sample = Dict{Any, Any}("t" => v)
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver,
		uneven_sampling = true,
		uneven_sampling_times = data_sample["t"])


	#println("Model: \n $model  \n ")

	#println("Data Sample: \n $data_sample, \n \n")
	#println("Measured Quantities: \n $measured_quantities \n ")

	for i in eachindex(measured_quantities)

		#        println(measured_quantities[i])
		#println(lhs(measured_quantities[i]))

		#        println(D(measured_quantities[i]))
		#        println(D(measured_quantities[i].lhs))
		#        println(D(measured_quantities[i].rhs))
		#        println(expand_derivatives(D(measured_quantities[i].lhs)))
		#        println(expand_derivatives(D(measured_quantities[i].rhs)))

		#       println(ModelingToolkit.diff2term(expand_derivatives(D(measured_quantities[i].lhs))))
		#       println(ModelingToolkit.diff2term(expand_derivatives(D(measured_quantities[i].rhs))))
		#    println(ModelingToolkit.diff2term(D(measured_quantities[i])))
	end

	sample_count = length(data_sample["t"])
	#   println(sample_count)

	#   println(model)
	#println("data_sample:")
	#println(data_sample)
	#println()


	SimpleParameterEstimation(model, measured_quantities, data_sample, solver)

end

main()

#res = ParameterEstimation.estimate(model, measured_quantities, data_sample, solver=solver)



