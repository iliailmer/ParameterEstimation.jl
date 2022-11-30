module ParameterEstimation
using Distributed
using ProgressMeter, Logging

using LinearAlgebra

using ModelingToolkit
import DifferentialEquations: Tsit5

using LinearSolve
using SIAN

using Nemo, Singular
# import Nemo: fmpq_mpoly, total_degree, vars, var_index, evaluate

using HomotopyContinuation

import TaylorSeries: Taylor1

using Groebner

using LinearAlgebra

Float = Union{Float64, Float32, Float16}

include("rational_interpolation/rational_interpolation.jl")
include("rational_interpolation/interpolant.jl")

include("identifiability/get_identifiability.jl")
include("identifiability/transcendence_basis.jl")
include("identifiability/identifiability_data.jl")
include("identifiability/utils.jl")

include("estimate.jl")
include("utils.jl")
include("metrics.jl")

export get_identifiability, estimate

end