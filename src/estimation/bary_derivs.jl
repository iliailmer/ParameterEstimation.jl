


#to use the below, you can just pass vectors of xvalues and yvalues like so:
#F = aaad(xdata, ydata)
#and then F will be a callable, i.e. F(0.5) should work.  Let's please restrict to real xdata, and if so it should be sorted.  F is only defined in the range of the xvalues.

#To construct a derivative, you can do
#derivf(z) = ForwardDiff.derivative(F, z)
#I hope and suspect that other AD frameworks should work as well.

function baryEval(z, f::Vector{T}, x::Vector{T}, w::Vector{T}, tol=1e-13) where {T}
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