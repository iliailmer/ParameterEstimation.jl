"""
	rational_interpolation_coefficients(x, y, n)
CODE COPIED FROM previous version of ParameterEstimation.jl
Perform a rational interpolation of the data `y` at the points `x` with numerator degree `n`.
This function only returns the coefficients of the numerator and denominator polynomials.

# Arguments
- `x`: the points where the data is sampled (e.g. time points).
- `y`: the data sample.
- `n`: the degree of the numerator.

# Returns
- `c`: the coefficients of the numerator polynomial.
- `d`: the coefficients of the denominator polynomial.
"""
function rational_interpolation_coefficients(x, y, n)
	N = length(x)
	m = N - n - 1
	A = zeros(N, N)
	if m > 0
		A_left_submatrix = reduce(hcat, [x .^ (i) for i in 0:(n)])
		A_right_submatrix = reduce(hcat, [x .^ (i) for i in 0:(m-1)])
		A = hcat(A_left_submatrix, -y .* A_right_submatrix)
		b = y .* (x .^ m)
		try
			prob = LinearSolve.LinearProblem(A, b)
			c = LinearSolve.solve(prob)
			return c[1:(n+1)], [c[(n+2):end]; 1]
		catch SingularException
			lu_res = lu(A)
			y = lu_res.L \ lu_res.P * b
			c = lu_res.U \ y
			return c[1:(n+1)], [c[(n+2):end]; 1]
		end

	else
		A = reduce(hcat, [x .^ i for i in 0:n])
		b = y
		prob = LinearSolve.LinearProblem(A, b)
		c = LinearSolve.solve(prob)
		return c, [1]
	end
end





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
		if ((z - x[j])^2 < sqrt(tol))
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

struct FHDapprox{T}
	internalFHD::T
end

(y::FHDapprox)(z) = baryEval(z, y.internalFHD.f, y.internalFHD.x, y.internalFHD.w)
(y::AAADapprox)(z) = baryEval(z, y.internalAAA.f, y.internalAAA.x, y.internalAAA.w)

function nth_deriv_at_deprecated(f, n::Int, t)  #todo(orebas) make this more efficient.
	if (n == 0)
		return f(t)
	elseif (n == 1)
		return ForwardDiff.derivative(f, t)
	else
		g(t) = nth_deriv_at(f, n - 1, t)
		return ForwardDiff.derivative(g, t)
	end
end


function nth_deriv_at(f, n::Int, t)  #todo(orebas) make this more efficient.
	if (n == 0)
		return f(t)
	else
		return TaylorDiff.derivative(f, t, n)
	end
end


function aaad(xs::AbstractArray{T}, ys::AbstractArray{T}) where {T}
	@suppress begin
		@assert length(xs) == length(ys)
		internalApprox = BaryRational.aaa(xs, ys)
		return AAADapprox(internalApprox)
	end
end



function fhd(xs::AbstractArray{T}, ys::AbstractArray{T}, N::Int) where {T}
	@suppress begin
		@assert length(xs) == length(ys)
		internalApprox = BaryRational.FHInterp(xs, ys, order = N)
		return FHDapprox(internalApprox)
	end
end

function fhdn(n)
	fh(xs, ys) = fhd(xs, ys, n)
	return fh
end


####################

struct FourierSeries
	m::Any
	b::Any
	K::Any
	cosines::Any
	sines::Any
end


function fourierEval(x, FS)
	z = FS.m * x + FS.b
	sum = 0.0
	for k in eachindex(FS.cosines)
		sum += FS.cosines[k] * cos((k) * z)
	end
	for k in eachindex(FS.sines)
		sum += FS.sines[k] * sin((k) * z)
	end
	sum += FS.K
	return sum
end

(y::FourierSeries)(z) = fourierEval(z, y)



function FourierInterp(xs, ys)
	@assert length(xs) == length(ys)
	N = length(xs)
	width = xs[end] - xs[begin]
	m = pi / width
	b = -pi * (xs[begin] / width + 0.5)
	f(t) = m * t + b
	scaledxs = f.(xs)
	sinescount = (N - 1) ÷ 2
	cosinescount = N - 1 - sinescount
	A = zeros(Float64, N, N)
	for i ∈ 1:N, j ∈ 1:N
		if (j == 1)
			A[i, 1] = 1
		elseif (j <= (cosinescount + 1))
			A[i, j] = cos((j - 1) * scaledxs[i])
		else
			temp = (j - cosinescount - 1)
			A[i, j] = sin(temp * scaledxs[i])
		end
	end

	prob = LinearProblem(A, ys, LinearSolve.OperatorCondition.VeryIllConditioned)
	sol = LinearSolve.solve(prob)
	X = sol.u
	temp = FourierSeries(m, b, X[begin], X[2:(cosinescount+1)], X[(cosinescount+2):end])
	return temp
end



struct BaryLagrange{T <: AbstractArray}
	x::T
	f::T
	w::T
end

function BarycentricLagrange(xs, ys)
	@assert length(xs) == length(ys)
	N = length(xs)
	w = ones(Float64, N)
	for k in eachindex(xs)
		for j in eachindex(xs)
			if (k != j)
				w[k] *= (xs[k] - xs[j])
			end
		end
		w[k] = 1 / w[k]
	end
	return BaryLagrange(xs, ys, w)
end

(y::BaryLagrange)(z) = baryEval(z, y.f, y.x, y.w)

struct RationalFunction{T <: AbstractArray}
	a::T
	b::T
end

(y::RationalFunction)(z) = rationaleval(z, y.a, y.b)
rationaleval(z, a, b) = evalpoly(z, a) / evalpoly(z, b)



function simpleratinterp(xs, ys, d1::Int)
	@assert length(xs) == length(ys)
	N = length(xs)
	d2 = N - d1 - 1
	A = zeros(Float64, N, N)
	for j in 1:N
		A[j, 1] = 1
		for k in 1:d1
			A[j, k+1] = xs[j]^k
		end
		for k in 1:d2
			A[j, d1+1+k] = -1.0 * ys[j] * xs[j]^k
		end
	end
	prob = LinearProblem(A, ys, LinearSolve.OperatorCondition.VeryIllConditioned)
	sol = LinearSolve.solve(prob)
	X = sol.u
	return RationalFunction(
		X[1:(d1+1)],
		[1; X[(d1+2):end]])
end






function betterratinterp(xs, ys, d1::Int)
	@assert length(xs) == length(ys)
	N = length(xs)
	(c, d) = rational_interpolation_coefficients(xs, ys, d1)
	return RationalFunction(c, d)
end



function SimpleRationalInterp(numerator_degree::Int)
	f(xs, ys) = simpleratinterp(xs, ys, numerator_degree)
	return f
end

function SimpleRationalInterpOld(numerator_degree::Int)
	f(xs, ys) = betterratinterp(xs, ys, numerator_degree)
	return f
end


function default_interpolator(datasize)
	interpolators = Dict(
			"AAA" => aaad,
			"FHD3" => fhdn(3),
			#			"Fourier" => FourierInterp,
			#			"BaryLagrange" => BarycentricLagrange,
		)
		if (datasize > 10)
			interpolators["FHD8"] = fhdn(8)
			#			interpolators["FHD6"] = fhdn(6)
			#stepsize = max(1, datasize ÷ 4)
			#for i in range(1, (datasize - 2), step = stepsize)
			#	interpolators["RatOld($i)"] = SimpleRationalInterpOld(i)
			#end
		end
return interpolators
end
