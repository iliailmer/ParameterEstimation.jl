using Test
using TestSetExtensions
using ModelingToolkit
using ParameterEstimation
using Nemo

@info "Testing started"

@testset "All the tests" begin @includetests ARGS end
