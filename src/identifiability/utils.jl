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

function count_solutions(identifiability_result::IdentifiabilityData)
    @info
    globally_id = identifiability_result.identifiability["globally"]
    locally_not_globally_id = identifiability_result.identifiability["locally_not_globally"]
    non_id = identifiability_result.identifiability["nonidentifiable"]
    weights_table = identifiability_result.weights
    solutions_table = Dict{Any, Int}()
    for param in globally_id
        solutions_table[param] = 1
    end
    basis_lex = Groebner.fglm(identifiability_result.basis)
    for param in locally_not_globally_id
    end
end