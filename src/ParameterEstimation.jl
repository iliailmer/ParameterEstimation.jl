module ParameterEstimation
import DifferentialEquations: Tsit5
import TaylorSeries: Taylor1

using ProgressMeter, Logging, Printf
using ModelingToolkit, LinearSolve, LinearAlgebra
using SIAN, Nemo, HomotopyContinuation, Groebner
using .ReturnCode

Float = Union{Float64, Float32, Float16}

include("rational_interpolation/rational_interpolation.jl")
include("rational_interpolation/interpolant.jl")

include("identifiability/check_identifiability.jl")
include("identifiability/transcendence_basis.jl")
include("identifiability/identifiability_data.jl")
include("identifiability/utils.jl")

include("result.jl")
include("estimate.jl")
include("filtering.jl")

include("utils.jl")
include("metrics.jl")

export check_identifiability, estimate, estimate_over_degrees, filter_solutions, sample_data

end