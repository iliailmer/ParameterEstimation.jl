"""
    check_identifiability(ode::ModelingToolkit.ODESystem; measured_quantities = Array{ModelingToolkit.Equation}[],
                          infolevel = 0)

Check identifiability of parameters in the ODE system `ode` using the
algorithm described in [1]. The function returns a `ParameterEstimation.IdentifiabilityData`
object that contains the results of the identifiability analysis.

# Arguments
    - `ode::ModelingToolkit.ODESystem`: The ODE system to be analyzed
    - `measured_quantities = Array{ModelingToolkit.Equation}[]`: A list of equations
        that define the measured quantities. If not provided, the outputs of the ODE
        system will be used.
    - `infolevel::Int`: The level of information to be printed during the analysis.

# References

[1] - Global Identifiability of Differential Models (Communications on Pure and Applied Mathematics, Volume 73, Issue 9, Pages 1831-1879, 2020.), https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.21921

[2] - SIAN: software for structural identifiability analysis of ODE models (Bioinformatics, Volume 35, Issue 16, Pages 2873--2874, 2019), https://academic.oup.com/bioinformatics/article/35/16/2873/5096881

[3] - https://github.com/alexeyovchinnikov/SIAN-Julia
"""
function check_identifiability(ode::ModelingToolkit.ODESystem;
                               measured_quantities = Array{ModelingToolkit.Equation}[],
                               infolevel = 0)
    if length(measured_quantities) == 0
        if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
            @info "Measured quantities are not provided, trying to find the outputs in input ODE."
            measured_quantities = filter(eq -> (ModelingToolkit.isoutput(eq.lhs)),
                                         ModelingToolkit.equations(ode))
        else
            throw(error("Measured quantities (output functions) were not provided and no outputs were found."))
        end
    end
    ode_prep, input_syms, gens_ = SIAN.PreprocessODE(ode, measured_quantities)
    t = ModelingToolkit.arguments(ModelingToolkit.states(ode)[1])[1]
    params_to_assess_ = SIAN.get_parameters(ode_prep)
    nemo2mtk = Dict(gens_ .=> input_syms)

    res = ParameterEstimation.identifiability_ode(ode_prep, params_to_assess_; p = 0.99,
                                                  p_mod = 0, infolevel = infolevel,
                                                  weighted_ordering = true,
                                                  local_only = false)

    @debug "Post-Processing: Converting Nemo output to ModelingToolkit types"
    out = Dict()
    for (id_type, pars) in pairs(res["identifiability"])
        out[id_type] = [ModelingToolkit.Num(substitute(nemo2mtk[each], t => 0))
                        for each in pars]
    end
    res["identifiability"] = out
    res["nemo_mtk"] = nemo2mtk
    return ParameterEstimation.IdentifiabilityData(res)
end

# Path: src/identifiability/identifiability_ode.jl
"""
    identifiability_ode(ode, params_to_assess; p = 0.99, p_mod = 0, infolevel = 0,
                        weighted_ordering = false, local_only = false)

Base function for identifiability assessment implmented in SIAN.jl. Detailed documentation is available in
`check_identifiability` function as well as in https://github.com/alexeyovchinnikov/SIAN-Julia.
"""
function identifiability_ode(ode, params_to_assess; p = 0.99, p_mod = 0, infolevel = 0,
                             weighted_ordering = false, local_only = false)
    @info "Solving the problem"

    if p_mod != 0
        @warn "Using `p_mod!=0` does not guarantee the same probability of correctness but allows to run the program faster. This warning was raised by `p_mod = $p_mod`."
    end
    # 1.Construct the maximal system

    # (a) ---------------

    @info "Constructing the maximal system"

    eqs, Q, x_eqs, y_eqs, x_vars, y_vars, u_vars, mu, all_indets, gens_Rjet = SIAN.get_equations(ode)

    non_jet_ring = ode.poly_ring
    z_aux = gens_Rjet[end - length(mu)]
    Rjet = gens_Rjet[1].parent

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    not_int_cond_params = gens_Rjet[(end - length(ode.parameters) + 1):end]
    all_params = vcat(not_int_cond_params, gens_Rjet[1:n])
    x_variables = gens_Rjet[1:n]
    for i in 1:(s + 1)
        x_variables = vcat(x_variables,
                           gens_Rjet[(i * (n + m + u) + 1):(i * (n + m + u) + n)])
    end
    u_variables = gens_Rjet[(n + m + 1):(n + m + u)]
    for i in 1:(s + 1)
        u_variables = vcat(u_variables,
                           gens_Rjet[((n + m + u) * i + n + m + 1):((n + m + u) * (i + 1))])
    end

    # (b,c) -------------
    X, X_eq = SIAN.get_x_eq(x_eqs, y_eqs, n, m, s, u, gens_Rjet)

    # (d,e) -----------
    Y, Y_eq = SIAN.get_y_eq(x_eqs, y_eqs, n, m, s, u, gens_Rjet)

    # 2.Truncate
    @info "Truncating"

    # (a) -----------------------
    d0 = BigInt(maximum(vcat([SIAN.Nemo.total_degree(SIAN.unpack_fraction(Q * eq[2])[1])
                              for eq in eqs], SIAN.Nemo.total_degree(Q))))

    # (b) -----------------------
    D1 = floor(BigInt,
               (length(params_to_assess) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) /
               (1 - p))

    # (c, d) ---------------
    sample = SIAN.sample_point(D1, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
    all_subs = sample[4]
    u_hat = sample[2]
    y_hat = sample[1]

    # (e) ------------------
    alpha = [1 for i in 1:n]
    beta = [0 for i in 1:m]
    Et = Array{SIAN.Nemo.fmpq_mpoly}(undef, 0)
    x_theta_vars = all_params
    prolongation_possible = [1 for i in 1:m]

    # (f) ------------------
    all_x_theta_vars_subs = SIAN.insert_zeros_to_vals(all_subs[1], all_subs[2])
    eqs_i_old = Array{SIAN.Nemo.fmpq_mpoly}(undef, 0)
    evl_old = Array{SIAN.Nemo.fmpq_mpoly}(undef, 0)
    while sum(prolongation_possible) > 0
        for i in 1:m
            if prolongation_possible[i] == 1
                eqs_i = vcat(Et, Y[i][beta[i] + 1])
                evl = [SIAN.Nemo.evaluate(eq, vcat(u_hat[1], y_hat[1]),
                                          vcat(u_hat[2], y_hat[2]))
                       for eq in eqs_i if !(eq in eqs_i_old)]
                evl_old = vcat(evl_old, evl)
                JacX = SIAN.jacobi_matrix(evl_old, x_theta_vars, all_x_theta_vars_subs)
                eqs_i_old = eqs_i
                if LinearAlgebra.rank(JacX) == length(eqs_i)
                    Et = vcat(Et, Y[i][beta[i] + 1])
                    beta[i] = beta[i] + 1
                    # adding necessary X-equations
                    polys_to_process = vcat(Et, [Y[k][beta[k] + 1] for k in 1:m])
                    while length(polys_to_process) != 0
                        new_to_process = Array{SIAN.Nemo.fmpq_mpoly}(undef, 0)
                        vrs = Set{SIAN.Nemo.fmpq_mpoly}()
                        for poly in polys_to_process
                            vrs = union(vrs,
                                        [v
                                         for v in SIAN.Nemo.vars(poly) if v in x_variables])
                        end
                        vars_to_add = Set{SIAN.Nemo.fmpq_mpoly}(v
                                                                for v in vrs
                                                                if !(v in x_theta_vars))
                        for v in vars_to_add
                            x_theta_vars = vcat(x_theta_vars, v)
                            ord_var = SIAN.get_order_var2(v, all_indets, n + m + u, s)
                            var_idx = SIAN.Nemo.var_index(ord_var[1])
                            poly = X[var_idx][ord_var[2]]
                            Et = vcat(Et, poly)
                            new_to_process = vcat(new_to_process, poly)
                            alpha[var_idx] = max(alpha[var_idx], ord_var[2] + 1)
                        end
                        polys_to_process = new_to_process
                    end
                else
                    prolongation_possible[i] = 0
                end
            end
        end
    end

    @info "Assessing local identifiability"

    max_rank = length(Et)
    for i in 1:m
        for j in (beta[i] + 1):length(Y[i])
            to_add = true
            for v in SIAN.get_vars(Y[i][j], x_vars, all_indets, n + m + u, s)
                if !(v in x_theta_vars)
                    to_add = false
                end
            end
            if to_add
                beta[i] = beta[i] + 1
                Et = vcat(Et, Y[i][j])
            end
        end
    end

    theta_l = Array{SIAN.Nemo.fmpq_mpoly}(undef, 0)
    params_to_assess_ = [SIAN.add_to_var(param, Rjet, 0) for param in params_to_assess]
    Et_eval_base = [SIAN.Nemo.evaluate(e, vcat(u_hat[1], y_hat[1]),
                                       vcat(u_hat[2], y_hat[2]))
                    for e in Et]
    for param_0 in params_to_assess_
        other_params = [v for v in x_theta_vars if v != param_0]
        Et_subs = [SIAN.Nemo.evaluate(e, [param_0],
                                      [SIAN.Nemo.evaluate(param_0, all_x_theta_vars_subs)])
                   for e in Et_eval_base]
        JacX = SIAN.jacobi_matrix(Et_subs, other_params, all_x_theta_vars_subs)
        if LinearAlgebra.rank(JacX) != max_rank
            theta_l = vcat(theta_l, param_0)
        end
    end
    x_theta_vars_reorder = vcat(theta_l,
                                reverse([x for x in x_theta_vars if !(x in theta_l)]))
    Et_ids, alg_indep = algebraic_independence(Et_eval_base, x_theta_vars_reorder,
                                               all_x_theta_vars_subs)
    @info "Found Pivots: [$(join(alg_indep, ", "))]"

    if length(theta_l) == 0
        @error "No identifiable parameters found!"
        @info "=== Summary ==="
        @info "Globally identifiable parameters:                 []"
        @info "Locally but not globally identifiable parameters: []"
        @info "Not identifiable parameters:                      [$(join(params_to_assess, ", "))]"

        return Dict("identifiability" => Dict("locally_identifiable" => [],
                                              "globally_identifiable" => [],
                                              "non_identifiable" => Set(SIAN.get_order_var(th,
                                                                                           non_jet_ring)[1]
                                                                        for th in params_to_assess)))
    else
        @info "Locally identifiable parameters: [$(join([SIAN.get_order_var(th, non_jet_ring)[1] for th in theta_l], ", "))]"
        @info "Not identifiable parameters:     [$(join([SIAN.get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess_, theta_l)], ", "))]"

        non_identifiable_parameters = Set(SIAN.get_order_var(th, non_jet_ring)[1]
                                          for th in setdiff(params_to_assess_, theta_l))

        if local_only
            return Dict("non_identifiable" => non_identifiable_parameters,
                        "locally_identifiable" => theta_l)
        end
        weights = Dict()
        if weighted_ordering
            weights = SIAN.get_weights(ode, non_identifiable_parameters)
        end
        # 3. Randomize.
        @info "Randomizing"
        # (a) ------------
        deg_variety = foldl(*, [BigInt(SIAN.Nemo.total_degree(e)) for e in Et])
        D2 = floor(BigInt,
                   3 / 4 * 6 * length(theta_l) * deg_variety *
                   (1 + 2 * d0 * maximum(beta)) /
                   (1 - p))
        # (b, c) ---------
        sample = SIAN.sample_point(D2, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq,
                                   Q)
        y_hat = sample[1]
        u_hat = sample[2]
        theta_hat = sample[3]

        # (d) ------------
        Et_hat = [SIAN.Nemo.evaluate(e, vcat(y_hat[1], u_hat[1]), vcat(y_hat[2], u_hat[2]))
                  for e in Et]
        transcendence_substitutions = Array{SIAN.Nemo.fmpq}(undef, 0)
        for (idx, var) in enumerate(theta_hat[1])
            if var in alg_indep
                push!(transcendence_substitutions, theta_hat[2][idx])
            end
        end
        Et_hat = [SIAN.Nemo.evaluate(e, alg_indep, transcendence_substitutions)
                  for e in Et_hat]
        Et_x_vars = Set{SIAN.Nemo.fmpq_mpoly}()
        for poly in Et_hat
            Et_x_vars = union(Et_x_vars, Set(SIAN.Nemo.vars(poly)))
        end
        Et_x_vars = setdiff(Et_x_vars, not_int_cond_params)
        Q_hat = SIAN.Nemo.evaluate(Q, u_hat[1], u_hat[2])
        vrs_sorted = vcat(sort([e for e in Et_x_vars],
                               lt = (x, y) -> SIAN.compare_diff_var(x, y, all_indets,
                                                                    n + m + u, s)), z_aux,
                          sort(not_int_cond_params, rev = true))

        # assign weights to variables
        if weighted_ordering
            for i in eachindex(Et_hat)
                for _var in Set(SIAN.Nemo.vars(Et_hat[i]))
                    _var_non_jet, _var_order = SIAN.get_order_var(_var, non_jet_ring)
                    Et_hat[i] = SIAN.make_substitution(Et_hat[i], _var,
                                                       _var^get(weights, _var_non_jet, 1),
                                                       parent(_var)(1))
                end
            end
            for _var in Set(SIAN.Nemo.vars(Q_hat))
                _var_non_jet, _var_order = SIAN.get_order_var(_var, non_jet_ring)
                Q_hat = SIAN.make_substitution(Q_hat, _var,
                                               _var^get(weights, _var_non_jet, 1),
                                               parent(_var)(1))
            end
        end
        # 4. Determine.
        @info "Gröbner basis computation"

        if p_mod > 0
            Et_hat = [SIAN._reduce_poly_mod_p(e, p_mod) for e in Et_hat]
            z_aux = SIAN._reduce_poly_mod_p(z_aux, p_mod)
            Q_hat = SIAN._reduce_poly_mod_p(Q_hat, p_mod)
            Rjet_new, vrs_sorted = SIAN.Nemo.PolynomialRing(SIAN.Nemo.GF(p_mod),
                                                            [string(v) for v in vrs_sorted],
                                                            ordering = :degrevlex)
        else
            Rjet_new, vrs_sorted = SIAN.Nemo.PolynomialRing(SIAN.Nemo.QQ,
                                                            [string(v) for v in vrs_sorted],
                                                            ordering = :degrevlex)
        end

        Et_hat = [SIAN.parent_ring_change(e, Rjet_new) for e in Et_hat]
        gb = groebner(vcat(Et_hat, SIAN.parent_ring_change(z_aux * Q_hat, Rjet_new) - 1))

        theta_g = Array{Any}(undef, 0)
        @info "Remainder computation"

        if p_mod > 0
            theta_l_new = [SIAN.parent_ring_change(SIAN._reduce_poly_mod_p(th, p_mod),
                                                   Rjet_new) for th in theta_l]

            for i in eachindex(theta_l)
                _var_non_jet, _var_order = SIAN.get_order_var(theta_l_new[i], non_jet_ring)
                if Groebner.normalform(gb, theta_l_new[i]^get(weights, _var_non_jet, 1)) ==
                   SIAN.parent_ring_change(SIAN._reduce_poly_mod_p(Rjet(theta_hat[2][findfirst(isequal(theta_l[i]),
                                                                                               theta_hat[1])]),
                                                                   p_mod), Rjet_new)
                    theta_g = vcat(theta_g, theta_l_new[i])
                end
            end
        else
            theta_l_new = [SIAN.parent_ring_change(th, Rjet_new) for th in theta_l]

            for i in eachindex(theta_l)
                _var_non_jet, _var_order = SIAN.get_order_var(theta_l_new[i], non_jet_ring)
                if Groebner.normalform(gb, theta_l_new[i]^get(weights, _var_non_jet, 1)) ==
                   SIAN.parent_ring_change(Rjet(theta_hat[2][findfirst(isequal(theta_l[i]),
                                                                       theta_hat[1])]),
                                           Rjet_new)
                    theta_g = vcat(theta_g, theta_l_new[i])
                end
            end
        end
        id_res_nemo = Dict("globally" => theta_g,
                           "locally_not_globally" => setdiff(theta_l_new, theta_g),
                           "nonidentifiable" => setdiff(params_to_assess_, theta_l))
        identifiability_result = Dict("globally" => Set(SIAN.get_order_var(th,
                                                                           non_jet_ring)[1]
                                                        for th in theta_g),
                                      "locally_not_globally" => Set(SIAN.get_order_var(th,
                                                                                       non_jet_ring)[1]
                                                                    for th in setdiff(theta_l_new,
                                                                                      theta_g)),
                                      "nonidentifiable" => Set(SIAN.get_order_var(th,
                                                                                  non_jet_ring)[1]
                                                               for th in setdiff(params_to_assess_,
                                                                                 theta_l)))
        @info "=== Summary ==="
        @info "Globally identifiable parameters:                 [$(join(identifiability_result["globally"], ", "))]"
        @info "Locally but not globally identifiable parameters: [$(join(identifiability_result["locally_not_globally"], ", "))]"
        @info "Not identifiable parameters:                      [$(join(identifiability_result["nonidentifiable"], ", "))]"
        @info "==============="

        # modify Y_eq for estimation purposes
        y_derivative_dict = Dict()
        for each in Y_eq
            name, order = SIAN.get_order_var(each[1], non_jet_ring)
            y_derivative_dict[each[1]] = order
        end
        Et_ = [Et[idx] for idx in Et_ids]
        full_result = Dict("full_polynomial_system" => [SIAN.Nemo.evaluate(e, alg_indep,
                                                                           transcendence_substitutions)
                                                        for e in Et],
                           "polynomial_system" => [SIAN.Nemo.evaluate(e, alg_indep,
                                                                      transcendence_substitutions)
                                                   for e in Et_],
                           "polynomial_system_to_solve" => HomotopyContinuation.System([]),
                           "denominator" => Q,
                           "Y_eq" => y_derivative_dict,
                           "vars" => vrs_sorted,
                           "vals" => all_subs,
                           "transcendence_basis_subs" => Dict(alg_indep .=>
                                                                  transcendence_substitutions),
                           "identifiability_nemo" => id_res_nemo,
                           "identifiability" => identifiability_result,
                           "basis" => gb,
                           "weights" => weights,
                           "non_jet_ring" => non_jet_ring)
        return full_result
    end
end
