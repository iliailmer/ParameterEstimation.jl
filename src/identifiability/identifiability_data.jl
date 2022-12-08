"""
    IdentifiabilityData

A struct that contains the data from identifiability analysis.
This is used for parameter estimation.
"""
struct IdentifiabilityData
    polynomial_system::Vector{Nemo.fmpq_mpoly}
    denomiantor::Nemo.fmpq_mpoly
    variables::Vector{Nemo.fmpq_mpoly}
    substitutions::Vector{Vector}
    identifiability_nemo::Any
    identifiability::Dict{Any, Any}
    transcendence_basis_subs::Vector{Nemo.RingElem}
    Y_eq::Dict{Any, Any}
    basis::Vector{Nemo.fmpq_mpoly}
    weights::Dict{fmpq_mpoly, Int64}
    non_jet_ring::Nemo.FmpqMPolyRing
    nemo_mtk::Dict
    solution_counts::Dict
    function IdentifiabilityData(input::Dict)
        solution_counts = count_solutions(input)
        return new(input["polynomial_system"], input["denominator"], input["vars"],
                   input["vals"], input["identifiability_nemo"], input["identifiability"],
                   input["transcendence_basis_subs"], input["Y_eq"], input["basis"],
                   input["weights"], input["non_jet_ring"], input["nemo_mtk"],
                   solution_counts)
    end
end
