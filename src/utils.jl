function nemo2hc(expr_tree::Union{Expr, Symbol})
    #traverse expr_tree
    if typeof(expr_tree) == Symbol
        return HomotopyContinuation.variables(expr_tree)[1]
    end
    if typeof(expr_tree) == Expr
        if expr_tree.head == :call
            if expr_tree.args[1] in [:+, :-, :*, :/, :^]
                return reduce(eval(expr_tree.args[1]), map(nemo2hc, expr_tree.args[2:end]))
            end
        end
    end
end

function nemo2hc(expr_tree::fmpq_mpoly)
    return nemo2hc(Meta.parse(string(expr_tree)))
end

function nemo2hc(expr_tree::Number)
    return expr_tree
end