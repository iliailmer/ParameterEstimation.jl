"""
    algebraic_independence(Et::Vector{SIAN.Nemo.fmpq_mpoly},
                           indets::Vector{SIAN.Nemo.fmpq_mpoly},
                           vals)

Returns the indices of the equations in Et to be used for polynomial solving
and the variables that form a transcendence basis.

# Arguments
- `Et::Vector{SIAN.Nemo.fmpq_mpoly}`: The equations to be solved (must come from identifiability check).
- `indets::Vector{SIAN.Nemo.fmpq_mpoly}`: The indeterminates.
- `vals::Vector{SIAN.Nemo.fmpq_mpoly}`: The values of the indeterminates sampled by identifiability algorithm.
"""
function algebraic_independence(Et::Vector{SIAN.Nemo.fmpq_mpoly},
                                indets::Vector{SIAN.Nemo.fmpq_mpoly},
                                vals)
    pivots = Vector{SIAN.Nemo.fmpq_mpoly}()
    Jacobian = SIAN.jacobi_matrix(Et, indets, vals)
    U = SIAN.Nemo.lu(Jacobian)[end]
    #find pivot columns in u
    for row_idx in 1:size(U, 1)
        row = U[row_idx, :]
        if !all(row .== 0)
            pivot_col = findfirst(row .!= 0)
            push!(pivots, indets[pivot_col[2]])
        end
    end
    current_idx = 1
    output_rows = Jacobian[current_idx, :]
    current_rank = 1
    output_ids = [1]
    for current_idx in 2:length(Et)
        current = [output_rows; Jacobian[current_idx, :]]
        if SIAN.Nemo.rank(current) > current_rank
            output_rows = current
            push!(output_ids, current_idx)
            current_rank += 1
        end
    end
    return output_ids, setdiff(indets, pivots)
end
