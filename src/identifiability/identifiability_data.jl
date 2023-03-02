PolySystem = Union{HomotopyContinuation.ModelKit.System, Vector{SIAN.Nemo.fmpq_mpoly}}

"""
    IdentifiabilityData

A struct that contains the data from identifiability analysis.
This is used for parameter estimation.

# Fields
- `polynomial_system::Vector{SIAN.Nemo.fmpq_mpoly}`: The polynomial system.
- `polynomial_system_to_solve::PolySystem`: The polynomial system with derivatives substitutited and ready to be solved.
- `denominator::SIAN.Nemo.fmpq_mpoly`: The denominator of the polynomial system.
- `variables::Vector{SIAN.Nemo.fmpq_mpoly}`: The variables of the polynomial system.
- `substitutions::Vector{Vector}`: The substitutions used to assess identifiability.
- `identifiability_nemo::Any`: The identifiability data from SIAN in Nemo data type.
- `identifiability::Dict`: The identifiability data from SIAN in HomotopyContinuation compatible data type.
- `basis::Vector{SIAN.Nemo.fmpq_mpoly}`: The transcendence basis of the polynomial system.
- `transcendence_basis_subs::Vector{SIAN.Nemo.RingElem}`: The transcendence basis substitutions of the polynomial system.
- `weights::Dict{SIAN.Nemo.fmpq_mpoly, Int64}`: The weights of the variables used by SIAN to assess GroebnerBasis.
"""
mutable struct IdentifiabilityData
    polynomial_system::Vector{SIAN.Nemo.fmpq_mpoly}
    polynomial_system_to_solve::PolySystem
    denomiantor::SIAN.Nemo.fmpq_mpoly
    variables::Vector{SIAN.Nemo.fmpq_mpoly}
    substitutions::Vector{Vector}
    identifiability_nemo::Any
    identifiability::AbstractDict
    transcendence_basis_subs::AbstractDict
    Y_eq::AbstractDict
    u_variables::AbstractDict
    basis::Vector{SIAN.Nemo.fmpq_mpoly}
    weights::AbstractDict{SIAN.Nemo.fmpq_mpoly, Int64}
    non_jet_ring::SIAN.Nemo.FmpqMPolyRing
    nemo_mtk::AbstractDict
    solution_counts::AbstractDict
    function IdentifiabilityData(input::AbstractDict)
        solution_counts = count_solutions(input)
        return new(input["polynomial_system"], input["polynomial_system_to_solve"],
                   input["denominator"], input["vars"],
                   input["vals"], input["identifiability_nemo"], input["identifiability"],
                   input["transcendence_basis_subs"], input["Y_eq"], input["u_variables"],
                   input["basis"],
                   input["weights"], input["non_jet_ring"], input["nemo_mtk"],
                   solution_counts)
    end
end
