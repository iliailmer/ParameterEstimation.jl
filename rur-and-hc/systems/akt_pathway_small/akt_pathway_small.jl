# Akt pathway
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L248C17-L279C10

# Akt pathway minus 8 parameters

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters reaction_1_k1 reaction_2_k2 reaction_3_k1 reaction_4_k1 reaction_5_k1 reaction_6_k1 reaction_7_k1 reaction_8_k1 reaction_9_k1
@variables t EGFR(t) pEGFR(t) pEGFR_Akt(t) Akt(t) pAkt(t) S6(t) pAkt_S6(t) pS6(t) EGF_EGFR(t) y1(t) y2(t) y3(t)
D = Differential(t)
states = [EGFR, pEGFR, pEGFR_Akt, Akt, pAkt, S6, pAkt_S6, pS6, EGF_EGFR]
parameters = [reaction_1_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1]

@info "Akt pathway: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(EGFR) ~
                0.1 * 0.1 + EGF_EGFR * 0.1 -
                EGFR * 0.1 - EGF_EGFR * reaction_1_k1,
            D(pEGFR) ~
                EGF_EGFR * reaction_9_k1 - pEGFR * reaction_4_k1 +
                pEGFR_Akt * reaction_2_k2 +
                pEGFR_Akt * reaction_3_k1 - Akt * pEGFR * 0.1,
            D(pEGFR_Akt) ~
                Akt * pEGFR * 0.1 - pEGFR_Akt * reaction_3_k1 -
                pEGFR_Akt * reaction_2_k2,
            D(Akt) ~
                pAkt * reaction_7_k1 + pEGFR_Akt * reaction_2_k2 -
                Akt * pEGFR * 0.1,
            D(pAkt) ~
                pAkt_S6 * 0.1 - pAkt * reaction_7_k1 +
                pAkt_S6 * reaction_6_k1 +
                pEGFR_Akt * reaction_3_k1 - S6 * pAkt * reaction_5_k1,
            D(S6) ~
                pAkt_S6 * 0.1 + pS6 * reaction_8_k1 -
                S6 * pAkt * reaction_5_k1,
            D(pAkt_S6) ~
                S6 * pAkt * reaction_5_k1 - pAkt_S6 * reaction_6_k1 -
                pAkt_S6 * 0.1,
            D(pS6) ~ pAkt_S6 * reaction_6_k1 - pS6 * reaction_8_k1,
            D(EGF_EGFR) ~
                EGF_EGFR * reaction_1_k1 - EGF_EGFR * reaction_9_k1 -
                EGF_EGFR * 0.1,
         ], t, states, parameters)
measured_quantities = [y1 ~ 0.1 * (pEGFR + pEGFR_Akt),
            y2 ~ 0.1 * (pAkt + pAkt_S6),
            y3 ~ 0.1 * pS6]

ic = [1.0 for i in 1:length(states)]
time_interval = [0.0, 1.0]
datasize = 21
p_true = [0.1 + i*(1 / (2length(parameters))) for i in 1:length(parameters)]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)

@info "" parameters p_true
@info "" data_sample

res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
        solver = solver)


