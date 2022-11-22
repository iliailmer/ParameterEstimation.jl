"""
    function algebraic_independence(Et::Vector{Nemo.fmpq_mpoly}, indets::Vector{Nemo.fmpq_mpoly})
"""
function algebraic_independence(Et::Vector{Nemo.fmpq_mpoly},
                                indets::Vector{Nemo.fmpq_mpoly},
                                vals)
    pivots = Vector{fmpq_mpoly}()
    Jacobian = SIAN.jacobi_matrix(Et, indets, vals)
    U = Nemo.lu(Jacobian)[end]
    #find pivot columns in u
    for row_idx in 1:size(U, 1)
        row = U[row_idx, :]
        if !all(row .== 0)
            pivot_col = findfirst(row .!= 0)
            push!(pivots, indets[pivot_col[2]])
        end
    end
    return setdiff(indets, pivots)
end