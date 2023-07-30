using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using BaryRational
using ForwardDiff



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

	#println("t")
	#println(t)

	#println("model_eq")
	#println(model_eq)

	#println("model_states")
	#println(model_states)

	#println("model_ps")
	#println(model_ps)

	t_vector = pop!(data_sample, "t") #TODO(orebas) make it use the independent variable name
	sample_count = length(t_vector)
	D = Differential(t)


	eqns = Array{Any}(undef, sample_count)

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
		#d2 = Dict()
		for j in eachindex(model_states)
			#dump(model_states[j])
			varname = (model_states[j].metadata.value[2])
			vname2 = Symbol("$(varname)_$(i)")
			#println(vname2)
			vname3 = @variables $vname2
			d[model_states[j]] = (vname3[1])

			dvar = ModelingToolkit.diff2term(expand_derivatives(D(model_states[j])))
			#println("DVAR ", dvar)
			dvarname = dvar.metadata.value[2]
			dvname2 = Symbol("$(dvarname)_$(i)")
			dvname3 = @variables $dvname2
			d[dvar] = dvname3[1]

		end

		#println(d)
		#for j in eachindex(model_states)
		#dump(model_states[j])
		#	varname = (model_states[j].metadata.value[2])
		#		vname2 = Symbol("$(varname)_$(i)")
		#			#println(vname2)
		#			d[model_states[j]] = (vname2)
		#			@variables $vname2
		#			dvar = ModelingToolkit.diff2term(expand_derivatives(D(model_states[j])))
		#			#println("DVAR ", dvar)
		#			dvarname = dvar.metadata.value[2]
		#			dvname2 = Symbol("$(dvarname)_$(i)")
		##			d[dvar] = dvname2

		#		end



		equations_time_slice_from_ODE_only = []# assemple the ODE equations
		#println(d)
		for j in eachindex(model_eq)
			println()
			println(model_eq[j])
			lhs1 = expand_derivatives(model_eq[j].lhs)
			rhs1 = expand_derivatives(model_eq[j].rhs)
			lhs2 = ModelingToolkit.diff2term(lhs1)
			rhs2 = ModelingToolkit.diff2term(rhs1)
			println()
			println(rhs2)
			println()
			println(d)
			lhs3 = substitute(lhs2, d)
			rhs3 = substitute(rhs2, d)
			#println(lhs1)
			#println(rhs1)
			#println(lhs2)
			#println(rhs2)
			println(lhs3)
			println(rhs3)
			push!(equations_time_slice_from_ODE_only, lhs3 ~ rhs3)
			#lhs = ModelingToolkit.diff2term(expand_derivatives(model_eq[j].lhs))
			#rhs = ModelingToolkit.diff2term(expand_derivatives(model_eq[j].rhs))
			#println("$i $j:")
			#println(lhs)
			#println(substitute(lhs, d))
			#println(rhs)
			#println(substitute(rhs, d))
		end
		equations_time_slice_from_measured_quantities_0th_deriv = []

		println()

		println()

		for j in measured_quantities
			r = j.rhs
			println("DEBUG1  ", r)
			r1 = ModelingToolkit.diff2term(expand_derivatives(r))
			println("DEBUG2  ", r1)
			#dump(r1)
			r2 = substitute(r1, d)
			#dump(r2)
			println("DEBUG3  ", r2)

			#yval = interpolants[r](data_sample[r][i])
			#eq = yval ~ r2
			#push!(equations_time_slice_from_measured_quantities_0th_deriv, eq)
		end
		println()

		println(equations_time_slice_from_ODE_only)
		println()

		println(equations_time_slice_from_measured_quantities_0th_deriv)
		println()


		#push!(eqns, model_converted)
		#println(measured_quantities)
		#println(data_sample)


	end
end

function main()
	solver = Tsit5()

	@parameters alpha
	@variables t theta_1(t) theta_1_dot(t) theta_2(t) theta_2_dot(t) y1(t) y2(t)
	D = Differential(t)

	#coupled pendulum equations
	states = [theta_1, theta_1_dot, theta_2, theta_2_dot]
	parameters = [alpha]
	@named model = ODESystem([
			D(theta_1) ~ theta_2 + theta_1,
			D(theta_1_dot) ~ -theta_1 * (alpha + 1) + alpha * theta_2,
			D(theta_2) ~ theta_2_dot,
			D(theta_2_dot) ~ alpha * theta_1 - theta_2 * (alpha + 1),
		], t, states, parameters)

	#initial conditions
	ic = [1.0, 0.0, 0.0, -1.0]
	p_true = [1.0]
	time_interval = [0.0, 10.0]
	datasize = 4

	v = randn(datasize)
	v = sort((v .- minimum(v)) / (maximum(v) - minimum(v))) * time_interval[2]

	measured_quantities = [y1 ~ theta_1, y2 ~ theta_2 + theta_1]
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

	a = theta_2
	b = alpha * theta_2
	c = theta_1 + theta_2
	d = a + b


	#println(a)
	#println(b)
	#println(c)
	#println(d)

	#println(typeof(a))
	#println(typeof(b))
	#println(typeof(c))
	#println(typeof(d))

	#dump(a)
	#dump(b)
	#dump(c)
	#dump(d)


	SimpleParameterEstimation(model, measured_quantities, data_sample, solver)

end

main()

#res = ParameterEstimation.estimate(model, measured_quantities, data_sample, solver=solver)



