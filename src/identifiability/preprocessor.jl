"""
    function preprocess_ode(de::ModelingToolkit.ODESystem, measured_quantities::Array{ModelingToolkit.Equation})

Input:
- `de` - ModelingToolkit.ODESystem, a system for identifiability query
- `measured_quantities` - array of output functions

Output:
- `ODE` object containing required data for identifiability assessment
"""
function preprocess_ode(de::ModelingToolkit.ODESystem,
                        measured_quantities::Array{ModelingToolkit.Equation},
                        inputs = Vector{Num}())
    @info "Preproccessing `ModelingToolkit.ODESystem` object"
    diff_eqs = filter(eq -> !(ModelingToolkit.isoutput(eq.lhs)),
                      ModelingToolkit.equations(de))
    y_functions = [each.lhs for each in measured_quantities]
    state_vars = filter(s -> !(ModelingToolkit.isinput(s) || ModelingToolkit.isoutput(s)),
                        ModelingToolkit.states(de))
    params = ModelingToolkit.parameters(de)
    t = ModelingToolkit.arguments(measured_quantities[1].lhs)[1]
    params_from_measured_quantities = ModelingToolkit.parameters(ModelingToolkit.ODESystem(measured_quantities,
                                                                                           t,
                                                                                           name = :DataSeries))
    params = union(params, params_from_measured_quantities)

    input_symbols = vcat(state_vars, y_functions, inputs, params)
    generators = string.(input_symbols)
    generators = map(g -> replace(g, "(t)" => ""), generators)
    R, gens_ = Nemo.PolynomialRing(Nemo.QQ, generators)
    state_eqn_dict = Dict{Oscar.fmpq_mpoly,
                          Union{Oscar.fmpq_mpoly,
                                Oscar.Generic.Frac{
                                                   Oscar.fmpq_mpoly
                                                   }}}()
    out_eqn_dict = Dict{Oscar.fmpq_mpoly,
                        Union{Oscar.fmpq_mpoly,
                              Oscar.Generic.Frac{
                                                 Oscar.fmpq_mpoly
                                                 }}}()

    for i in eachindex(diff_eqs)
        if !(typeof(diff_eqs[i].rhs) <: Number)
            state_eqn_dict[substitute(state_vars[i], input_symbols .=> gens_)] = eval_at_nemo(diff_eqs[i].rhs,
                                                                                              Dict(input_symbols .=>
                                                                                                       gens_))
        else
            state_eqn_dict[substitute(state_vars[i], input_symbols .=> gens_)] = R(diff_eqs[i].rhs)
        end
    end
    for i in 1:length(measured_quantities)
        out_eqn_dict[substitute(y_functions[i], input_symbols .=> gens_)] = eval_at_nemo(measured_quantities[i].rhs,
                                                                                         Dict(input_symbols .=>
                                                                                                  gens_))
    end

    inputs_ = [substitute(each, input_symbols .=> gens_) for each in inputs]
    if isequal(length(inputs_), 0)
        inputs_ = Vector{Oscar.fmpq_mpoly}()
    end
    return (ODE{Oscar.fmpq_mpoly}(state_eqn_dict,
                                  out_eqn_dict,
                                  inputs_),
            input_symbols, gens_)
end
