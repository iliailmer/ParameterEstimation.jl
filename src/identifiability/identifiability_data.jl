struct IdentifiabilityData
    polynomial_system::Vector{Nemo.fmpq_mpoly}
    denomiantor::Nemo.fmpq_mpoly
    variables::Vector{Nemo.fmpq_mpoly}
    substitutions::Vector{Vector}
    identifiability::Dict{Any, Any}
    transcendence_basis_subs::Vector{Nemo.RingElem}
    Y_eq::Dict{Any, Any}
    basis::Vector{Nemo.fmpq_mpoly}
    weights::Dict{fmpq_mpoly, Int64}
    non_jet_ring::Nemo.FmpqMPolyRing
    function IdentifiabilityData(input::Dict)
        return new(input["polynomial_system"], input["denominator"], input["vars"],
                   input["vals"],
                   input["identifiability"], input["transcendence_basis_subs"],
                   input["Y_eq"], input["basis"], input["weights"], input["non_jet_ring"])
    end
end
