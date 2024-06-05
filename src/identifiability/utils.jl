function Base.show(io::IO, identifiability_result::IdentifiabilityData)
    println("=== Summary ===")
    println("Globally identifiable parameters:                 [$(join(identifiability_result.identifiability["globally"], ", "))]")
    println("Locally but not globally identifiable parameters: [$(join(identifiability_result.identifiability["locally_not_globally"], ", "))]")
    println("Not identifiable parameters:                      [$(join(identifiability_result.identifiability["nonidentifiable"], ", "))]")
    println("===============")
end

function Base.getindex(identifiability_result::IdentifiabilityData, key::Symbol)
    return getproperty(identifiability_result, key)
end

function Base.getindex(identifiability_result::IdentifiabilityData, key::String)
    return getproperty(identifiability_result, Symbol(key))
end

function Base.setindex!(identifiability_result::IdentifiabilityData,
        value::PolySystem, key::String)
    return setproperty!(identifiability_result, Symbol(key), value)
end

function Base.setindex!(identifiability_result::IdentifiabilityData,
        value::PolySystem, key::Symbol)
    return setproperty!(identifiability_result, key, value)
end

"""
    count_solutions(identifiability_result::IdentifiabilityData)

Count the number of solutions for each parameter. This uses identifiability result from `identifiability_analysis`.
The result corresponds to number of solutions for each locally identifiable parameter.

Non-identifable parameters are assigned 0. Globally identifable parameters are assigned 1.

# Arguments
    - `identifiability_result::IdentifiabilityData`: Identifiability result from `identifiability_analysis`.
"""
function count_solutions(identifiability_result)
    throw("Not implemented, sorry.")
    # globally_id = identifiability_result["identifiability_nemo"]["globally"]
    # locally_not_globally_id = identifiability_result["identifiability_nemo"]["locally_not_globally"]
    # non_id = identifiability_result["identifiability_nemo"]["nonidentifiable"]
    # weights_table = identifiability_result["weights"]
    # non_jet_ring = identifiability_result["non_jet_ring"]
    # n2m = identifiability_result["nemo_mtk"]
    # solutions_table = Dict{Any, Int}()
    # for param in globally_id
    #     v = SIAN.get_order_var(param, non_jet_ring)[1]
    #     solutions_table[n2m[v]] = 1
    # end
    # for param in non_id
    #     v = SIAN.get_order_var(param, non_jet_ring)[1]
    #     solutions_table[n2m[v]] = 0
    # end
    # # basis_lex = Groebner.fglm(identifiability_result.basis)
    # R = parent(identifiability_result["basis"][1])
    # OR, svars = Oscar.PolynomialRing(Oscar.QQ, string.(gens(R)); ordering = :lex)
    # basis_sing = [SIAN.parent_ring_change(identifiability_result["basis"][i], OR)
    #               for i in 1:length(identifiability_result["basis"])]
    # locally_not_globally_id_sing = [SIAN.parent_ring_change(locally_not_globally_id[i],
    #                                                         OR)
    #                                 for i in eachindex(locally_not_globally_id)]
    # OIdeal = Oscar.ideal(OR, basis_sing)

    # for (idx, param) in enumerate(locally_not_globally_id_sing)
    #     others = setdiff(svars, [param])
    #     elim = Oscar.gens(Oscar.eliminate(OIdeal, others))
    #     polynomials = filter(x -> param in Oscar.vars(x),
    #                          elim)[1]
    #     param_idx = findfirst(x -> x == param, svars)
    #     v = SIAN.get_order_var(locally_not_globally_id[idx], non_jet_ring)[1]
    #     dgrs = Oscar.degrees(polynomials)[param_idx] /
    #            get(weights_table, v, 1)
    #     solutions_table[n2m[v]] = Int(dgrs)
    # end
    # return solutions_table
end
