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
    state_vars = ModelingToolkit.unknowns(de)
    params = ModelingToolkit.parameters(de)
    t = ModelingToolkit.arguments(measured_quantities[1].lhs)[1]
    params_from_measured_quantities = ModelingToolkit.parameters(ModelingToolkit.ODESystem(
        measured_quantities,
        t,
        name = :DataSeries))
    params = union(params, params_from_measured_quantities)

    input_symbols = vcat(state_vars, y_functions, inputs, params)
    generators = string.(input_symbols)
    generators = map(g -> replace(g, "(t)" => ""), generators)
    R, gens_ = Nemo.polynomial_ring(Nemo.QQ, generators)
    state_eqn_dict = Dict{Nemo.QQMPolyRingElem,
        Union{Nemo.QQMPolyRingElem,
            Nemo.Generic.FracFieldElem{
                Nemo.QQMPolyRingElem
            }}}()
    out_eqn_dict = Dict{Nemo.QQMPolyRingElem,
        Union{Nemo.QQMPolyRingElem,
            Nemo.Generic.FracFieldElem{
                Nemo.QQMPolyRingElem
            }}}()

    for i in eachindex(diff_eqs)
        if !(typeof(diff_eqs[i].rhs) <: Number)
            state_eqn_dict[substitute(state_vars[i], input_symbols .=> gens_)] = StructuralIdentifiability.eval_at_nemo(
                diff_eqs[i].rhs,
                Dict(input_symbols .=>
                    gens_))
        else
            state_eqn_dict[substitute(state_vars[i], input_symbols .=> gens_)] = R(diff_eqs[i].rhs)
        end
    end
    for i in 1:length(measured_quantities)
        out_eqn_dict[substitute(y_functions[i], input_symbols .=> gens_)] = StructuralIdentifiability.eval_at_nemo(
            measured_quantities[i].rhs,
            Dict(input_symbols .=>
                gens_))
    end

    inputs_ = [substitute(each, input_symbols .=> gens_) for each in inputs]
    if isequal(length(inputs_), 0)
        inputs_ = Vector{Nemo.QQMPolyRingElem}()
    end
    return (ODE{Nemo.QQMPolyRingElem}(state_eqn_dict,
            out_eqn_dict,
            inputs_),
        input_symbols, gens_)
end
