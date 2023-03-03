var documenterSearchIndex = {"docs":
[{"location":"library/identifiability/identifiability/#Identifiability-Information","page":"Identifiability","title":"Identifiability Information","text":"","category":"section"},{"location":"library/identifiability/identifiability/","page":"Identifiability","title":"Identifiability","text":"ParameterEstimation.check_identifiability","category":"page"},{"location":"library/identifiability/identifiability/#ParameterEstimation.check_identifiability","page":"Identifiability","title":"ParameterEstimation.check_identifiability","text":"function checkidentifiability(ode::ModelingToolkit.ODESystem;                                measuredquantities = Array{ModelingToolkit.Equation}[],                                inputs::Vector{Num} = Array{Num}[],                                infolevel = 0)\n\nCheck identifiability of parameters in the ODE system ode using the algorithm described in [1]. The function returns a ParameterEstimation.IdentifiabilityData object that contains the results of the identifiability analysis.\n\nArguments\n\n- `ode::ModelingToolkit.ODESystem`: The ODE system to be analyzed\n- `measured_quantities = Array{ModelingToolkit.Equation}[]`: A list of equations\n    that define the measured quantities. If not provided, the outputs of the ODE\n    system will be used.\n- `inputs::Vector{Num} = Array{Num}[]`: A list of input functions, if any are present.\n- `infolevel::Int`: The level of information to be printed during the analysis.\n\nReferences\n\n[1] - Global Identifiability of Differential Models (Communications on Pure and Applied Mathematics, Volume 73, Issue 9, Pages 1831-1879, 2020.), https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.21921\n\n[2] - SIAN: software for structural identifiability analysis of ODE models (Bioinformatics, Volume 35, Issue 16, Pages 2873–2874, 2019), https://academic.oup.com/bioinformatics/article/35/16/2873/5096881\n\n[3] - https://github.com/alexeyovchinnikov/SIAN-Julia\n\n\n\n\n\n","category":"function"},{"location":"library/identifiability/identifiability/","page":"Identifiability","title":"Identifiability","text":"ParameterEstimation.IdentifiabilityData","category":"page"},{"location":"library/identifiability/identifiability/#ParameterEstimation.IdentifiabilityData","page":"Identifiability","title":"ParameterEstimation.IdentifiabilityData","text":"IdentifiabilityData\n\nA struct that contains the data from identifiability analysis. This is used for parameter estimation.\n\nFields\n\npolynomial_system::Vector{SIAN.Nemo.fmpq_mpoly}: The polynomial system.\npolynomial_system_to_solve::PolySystem: The polynomial system with derivatives substitutited and ready to be solved.\ndenominator::SIAN.Nemo.fmpq_mpoly: The denominator of the polynomial system.\nvariables::Vector{SIAN.Nemo.fmpq_mpoly}: The variables of the polynomial system.\nsubstitutions::Vector{Vector}: The substitutions used to assess identifiability.\nidentifiability_nemo::Any: The identifiability data from SIAN in Nemo data type.\nidentifiability::Dict: The identifiability data from SIAN in HomotopyContinuation compatible data type.\nbasis::Vector{SIAN.Nemo.fmpq_mpoly}: The transcendence basis of the polynomial system.\ntranscendence_basis_subs::Vector{SIAN.Nemo.RingElem}: The transcendence basis substitutions of the polynomial system.\nweights::Dict{SIAN.Nemo.fmpq_mpoly, Int64}: The weights of the variables used by SIAN to assess GroebnerBasis.\n\n\n\n\n\n","category":"type"},{"location":"tutorials/estimate/#Parameter-Estimation","page":"Simple Example","title":"Parameter Estimation","text":"","category":"section"},{"location":"tutorials/estimate/#Introduction","page":"Simple Example","title":"Introduction","text":"","category":"section"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"In this tutorial, we provide a general overview of using ParameterEstimation.jl.","category":"page"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"Assume we have a simple ODE model with output as below","category":"page"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"begincasesdotx = -mu xy = x^2+xendcases","category":"page"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"If we collect the sample at 4 time points between 0 and 1, we obtain a collection:","category":"page"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"\"t\"     => [0.000, 0.333, 0.666, 1.000],\n  x^2 + x => [2.000, 1.563, 1.229, 0.974]","category":"page"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"This is all that is needed for the program: a symbolic model (ODE and outputs) and a dictionary of data.","category":"page"},{"location":"tutorials/estimate/#Code","page":"Simple Example","title":"Code","text":"","category":"section"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"Below is the working code example:","category":"page"},{"location":"tutorials/estimate/","page":"Simple Example","title":"Simple Example","text":"using ParameterEstimation\nusing ModelingToolkit\n\n# Input:\n# -- Differential model\n@parameters mu\n@variables t x(t) y(t)\nD = Differential(t)\n@named Sigma = ODESystem([D(x) ~ -mu * x],\n                         t, [x], [mu])\nouts = [y ~ x^2 + x]\n\n# -- Data\ndata = Dict(\n  \"t\"     => [0.000, 0.333, 0.666, 1.000],\n  x^2 + x => [2.000, 1.563, 1.229, 0.974])\n\n# Run\nres = estimate(Sigma, outs, data);","category":"page"},{"location":"library/rational_interpolation/interpolant/#Interpolant","page":"Interpolant","title":"Interpolant","text":"","category":"section"},{"location":"library/rational_interpolation/interpolant/","page":"Interpolant","title":"Interpolant","text":"ParameterEstimation.Interpolant","category":"page"},{"location":"library/rational_interpolation/interpolant/#ParameterEstimation.Interpolant","page":"Interpolant","title":"ParameterEstimation.Interpolant","text":"Interpolant\n\nA structure that stores information about the interpolation result.\n\nFields\n\nI: the single variable differentiable function, result of interpolation.\ndIdt: the TaylorSeries derivative of I.\n\n\n\n\n\n","category":"type"},{"location":"library/rational_interpolation/rational_interpolation/#Rational-Interpolation-utils","page":"Rational Interpolation","title":"Rational Interpolation utils","text":"","category":"section"},{"location":"library/rational_interpolation/rational_interpolation/","page":"Rational Interpolation","title":"Rational Interpolation","text":"ParameterEstimation.interpolate","category":"page"},{"location":"library/rational_interpolation/rational_interpolation/#ParameterEstimation.interpolate","page":"Rational Interpolation","title":"ParameterEstimation.interpolate","text":"interpolate(identifiability_result, data_sample,\n            measured_quantities; interpolation_degree::Int = 1,\n            diff_order::Int = 1, at_t::Float = 0.0,\n            method::Symbol = :homotopy)\n\nThis function performs the key step in parameter estimation.\n\nIt interpolates the data in `data_sample` and computes the `TaylorSeries` expansion.\nThese results are stored in the `Interpolant` object and are applied to the polynomial system in `identifiability_result`.\n\nArguments\n\nidentifiability_result: the result of the identifiability check.\ndata_sample: a dictionary of the data samples. The keys are the symbols of the measured quantities and the values are the data samples.\nmeasured_quantities: the measured quantities (outputs as equations of the form y ~ x).\ninterpolation_degree::Int = 1: the degree of the numerator of the rational interpolation.\ndiff_order::Int = 1: the order of the derivative to be computed.\nat_t::Float = 0.0: the time point where the Taylor series expansion is computed.\nmethod::Symbol = :homotopy: the method used to solve the polynomial system. Can be one of :homotopy (recommended) and :msolve.\n\nReturns\n\nSystem: the polynomial system with the interpolated data applied. This system is compatible with HomotopyContinuation solving.\n\n\n\n\n\ninterpolate(time, sample, numer_degree::Int, diff_order::Int = 1, at_t::Float = 0.0)\n\nThis function performs a rational interpolation of the data sample at the points time with numerator degree numer_degree. It returns an Interpolant object that contains the interpolated function and its derivatives.\n\n\n\n\n\n","category":"function"},{"location":"library/rational_interpolation/rational_interpolation/","page":"Rational Interpolation","title":"Rational Interpolation","text":"ParameterEstimation.rational_interpolation_coefficients","category":"page"},{"location":"library/rational_interpolation/rational_interpolation/#ParameterEstimation.rational_interpolation_coefficients","page":"Rational Interpolation","title":"ParameterEstimation.rational_interpolation_coefficients","text":"rational_interpolation_coefficients(x, y, n)\n\nPerform a rational interpolation of the data y at the points x with numerator degree n. This function only returns the coefficients of the numerator and denominator polynomials.\n\nArguments\n\nx: the points where the data is sampled (e.g. time points).\ny: the data sample.\nn: the degree of the numerator.\n\nReturns\n\nc: the coefficients of the numerator polynomial.\nd: the coefficients of the denominator polynomial.\n\n\n\n\n\n","category":"function"},{"location":"library/estimate/#Estimate","page":"Estimation","title":"Estimate","text":"","category":"section"},{"location":"library/estimate/","page":"Estimation","title":"Estimation","text":"ParameterEstimation.estimate","category":"page"},{"location":"library/estimate/#ParameterEstimation.estimate","page":"Estimation","title":"ParameterEstimation.estimate","text":"estimate(model::ModelingToolkit.ODESystem,\n        measured_quantities::Vector{ModelingToolkit.Equation},\n        data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}();\n        at_time::T = 0.0, method = :homotopy, solver = Tsit5(),\n        degree_range = nothing, real_tol = 1e-10,\n        threaded = Threads.nthreads() > 1) where {T <: Float}\n\nRun estimation over a range of interpolation degrees. Return the best estimate according to a heuristic:     - the best estimate is the one with the smallest error between sample data and ODE solution with current parameters (estimates);\n\nArguments\n\nmodel::ModelingToolkit.ODESystem: the ODE model;\nmeasured_quantities::Vector{ModelingToolkit.Equation}: the measured quantities (output functions that were sampled experimentally);\ndata_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(): the data sample, a dictionary with keys being the measured quantities and                                                               values being the corresponding data. Must include the time vector;\nat_time::T = 0.0: the time used for derivative computation;\nmethod = :homotopy: the method used for polynomial system solving. Can be one of :homotopy (recommended) or :msolve;\nsolver: the ODE solver used for ODE solution computation (default: Tsit5());\ndegree_range = nothing: the range of interpolation degrees to be used. If nothing, the range is computed automatically;\nreal_tol = 1e-10: the tolerance used for real root finding;\nthreaded = Threads.nthreads() > 1: whether to use multiple threads for computation (determined automatically).\n\nReturns\n\nresult::Vector{EstimationResult}: the result of the estimation, a vector of EstimationResult objects.\n\n\n\n\n\n","category":"function"},{"location":"library/estimate/","page":"Estimation","title":"Estimation","text":"ParameterEstimation.estimate","category":"page"},{"location":"library/estimate/","page":"Estimation","title":"Estimation","text":"ParameterEstimation.EstimationResult","category":"page"},{"location":"library/estimate/#ParameterEstimation.EstimationResult","page":"Estimation","title":"ParameterEstimation.EstimationResult","text":"EstimationResult\n\nA container for the results of an estimation. Contains the estimated parameters and initial conditions (state values at a given time), the degree of the rational interpolation used, the error between the estimated ODE solution and the sample data, and the return code.\n\nFields\n\nparameters::OrderedDict: The estimated parameters.\nstates::OrderedDict: The estimated initial conditions.\ndegree::Int64: The degree of the rational interpolation used.\nat_time::Float64: The time at which the initial conditions are estimated.\nerr::Union{Nothing, Float64}: The error between the estimated ODE solution and the sample data.\ninterpolants::Union{Nothing, Dict{Any, Interpolant}}: The rational interpolants used to estimate the parameters and initial conditions.\nreturn_code::Any: The return code of the estimation.\ndatasize::Int64: The number of data points used in the estimation.\n\n\n\n\n\n","category":"type"},{"location":"#ParameterEstimation.jl","page":"Home","title":"ParameterEstimation.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Software for Parameter Estimation Based on Identifiability Information and Real Data","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install ParameterEstimation.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"ParameterEstimation\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"To","category":"page"},{"location":"#Citation","page":"Home","title":"Citation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TBA","category":"page"},{"location":"#Feature-Summary","page":"Home","title":"Feature Summary","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Parameter estimation based on sample data\nEstimated values are reported based on identifiability: local (finitely many), global (unque), unidentifiable.","category":"page"},{"location":"library/filtering/#Filter-Estimated-Values","page":"Filtering","title":"Filter Estimated Values","text":"","category":"section"},{"location":"library/filtering/","page":"Filtering","title":"Filtering","text":"ParameterEstimation.filter_solutions","category":"page"},{"location":"library/filtering/#ParameterEstimation.filter_solutions","page":"Filtering","title":"ParameterEstimation.filter_solutions","text":"filter_solutions(results::Vector{EstimationResult},\n                identifiability_result::IdentifiabilityData,\n                model::ModelingToolkit.ODESystem,\n                inputs::Vector{ModelingToolkit.Equation},\n                data_sample::AbstractDict{Any, Vector{T}} = Dict{Any, Vector{T}}();\n                solver = Tsit5(),\n                topk = 1) where {T <: Float}\n\nFilter estimation results stored in results vector based on ODE solving and checking against the sample. In addition, takes into account global and local identifiability of parameters when filtering.\n\nArguments\n\nresults::Vector{EstimationResult}: the vector of estimation results.\nidentifiability_result::IdentifiabilityData: the result of identifiability analysis.\nmodel::ModelingToolkit.ODESystem: the ODE system.\ninputs::Vector{ModelingToolkit.Equation}: the inputs of the ODE system.\ndata_sample::AbstractDict{Any, Vector{T}} = Dict{Any, Vector{T}}(): the data sample used for estimation (same functions as measured_quantities).                                                               The keys of the dictionary are the measured quantities                                                               and the values are the corresponding data samples.\ntime_interval::Vector{T} = Vector{T}(): the time interval of the ODE system.\ntopk = 1: (optional) the number of best estimates to return.\n\nReturns\n\nEstimationResult: the best estimate (if topk = 1) or the vector of best estimates (if topk > 1).\n\n\n\n\n\n","category":"function"},{"location":"library/filtering/","page":"Filtering","title":"Filtering","text":"ParameterEstimation.solve_ode\nParameterEstimation.solve_ode!","category":"page"},{"location":"library/filtering/#ParameterEstimation.solve_ode","page":"Filtering","title":"ParameterEstimation.solve_ode","text":"solve_ode(model, estimate::EstimationResult, data_sample; solver = Tsit5(),\n          return_ode = false)\n\nSolves the ODE system model with the parameters and initial conditions given by estimate. Compute the error between the solution and the data sample. The error is recorded in the EstimationResult.\n\nArguments\n\nmodel: the ODE system to be solved.\nestimate::EstimationResult: the parameters and initial conditions of the ODE system.\ndata_sample: the data sample used for estimation (same functions as measured_quantities).                The keys of the dictionary are the measured quantities                and the values are the corresponding data samples.\nsolver = Tsit5(): (optional) the solver used to solve the ODE system, see DifferentialEquations for available solvers.\nreturn_ode = false: (optional) whether to return the ODE solution.\n\nReturns\n\nodesolution: the solution of the ODE system (if `returnodeis set totrue`).\nEstimationResult: the estimated parameters and initial conditions of the model.\n\n\n\n\n\n","category":"function"},{"location":"library/filtering/#ParameterEstimation.solve_ode!","page":"Filtering","title":"ParameterEstimation.solve_ode!","text":"solve_ode!(model, estimates::Vector{EstimationResult},  inputs::Vector{Equation},data_sample; solver = Tsit5())\n\nRun solve_ode for multiple estimates and store the results (error between solution and sample) in each estimate. This is done in-place.\n\n\n\n\n\n","category":"function"}]
}
