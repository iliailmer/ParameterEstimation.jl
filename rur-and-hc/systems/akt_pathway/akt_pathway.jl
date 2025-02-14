# Akt pathway
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L248C17-L279C10

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters pro_EGFR EGFR_turnover a1 a2 a3 reaction_1_k1 reaction_1_k2 reaction_2_k1 reaction_2_k2 reaction_3_k1 reaction_4_k1 reaction_5_k1 reaction_5_k2 reaction_6_k1 reaction_7_k1 reaction_8_k1 reaction_9_k1
@variables t EGFR(t) pEGFR(t) pEGFR_Akt(t) Akt(t) pAkt(t) S6(t) pAkt_S6(t) pS6(t) EGF_EGFR(t) y1(t) y2(t) y3(t)
D = Differential(t)
states = [EGFR, pEGFR, pEGFR_Akt, Akt, pAkt, S6, pAkt_S6, pS6, EGF_EGFR]
parameters = [pro_EGFR, EGFR_turnover, a1, a2, a3, reaction_1_k1, reaction_1_k2, reaction_2_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_5_k2, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1]

@info "Akt pathway: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(EGFR) ~
                EGFR_turnover * pro_EGFR + EGF_EGFR * reaction_1_k2 -
                EGFR * EGFR_turnover - EGF_EGFR * reaction_1_k1,
            D(pEGFR) ~
                EGF_EGFR * reaction_9_k1 - pEGFR * reaction_4_k1 +
                pEGFR_Akt * reaction_2_k2 +
                pEGFR_Akt * reaction_3_k1 - Akt * pEGFR * reaction_2_k1,
            D(pEGFR_Akt) ~
                Akt * pEGFR * reaction_2_k1 - pEGFR_Akt * reaction_3_k1 -
                pEGFR_Akt * reaction_2_k2,
            D(Akt) ~
                pAkt * reaction_7_k1 + pEGFR_Akt * reaction_2_k2 -
                Akt * pEGFR * reaction_2_k1,
            D(pAkt) ~
                pAkt_S6 * reaction_5_k2 - pAkt * reaction_7_k1 +
                pAkt_S6 * reaction_6_k1 +
                pEGFR_Akt * reaction_3_k1 - S6 * pAkt * reaction_5_k1,
            D(S6) ~
                pAkt_S6 * reaction_5_k2 + pS6 * reaction_8_k1 -
                S6 * pAkt * reaction_5_k1,
            D(pAkt_S6) ~
                S6 * pAkt * reaction_5_k1 - pAkt_S6 * reaction_6_k1 -
                pAkt_S6 * reaction_5_k2,
            D(pS6) ~ pAkt_S6 * reaction_6_k1 - pS6 * reaction_8_k1,
            D(EGF_EGFR) ~
                EGF_EGFR * reaction_1_k1 - EGF_EGFR * reaction_9_k1 -
                EGF_EGFR * reaction_1_k2,
         ], t, states, parameters)
measured_quantities = [y1 ~ a1 * (pEGFR + pEGFR_Akt),
            y2 ~ a2 * (pAkt + pAkt_S6),
            y3 ~ a3 * pS6]

ic = [1.0 for i in 1:length(states)]
time_interval = [0.0, 1.0]
datasize = 21
p_true = [0.1 + i*(1 / (2length(parameters))) for i in 1:length(parameters)]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)

@info "" data_sample
println(parameters)
println(p_true)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
        solver = solver)


