amigo_dir = "./amigo/"
param_estim_dir = "./pe.jl/"
sciml_dir = "./sciml/"
iqm_dir = "./iqm/"
using Plots
using DataFrames
using CSV

function read_csv_data(filename)
    data = CSV.read(filename, DataFrame, delim = ",")
    return data
end

function read_txt_data(filename)
    data = CSV.read(filename, DataFrame, delim = " ",
                    header = ["num_points", "max_abs_rel_err"])
    return data
end

amigo = readdir(amigo_dir)
pe = readdir(param_estim_dir)
sciml = readdir(sciml_dir)
iqm = readdir(iqm_dir)

model_name_dict = Dict("FHN_error.csv" => ["Fitzhugh-Nagumo Model", 5],
                       "HIV_error.csv" => ["HIV Model", 15],
                       "crauste_error.csv" => ["Crauste Model", 18],
                       "daisy_ex3_error.csv" => ["DAISY 3-Compartment Model 1", 9],
                       "daisy_mamil3_error.csv" => ["DAISY 3-Compartment Model 2", 8],
                       "daisy_mamil4_error.csv" => ["DAISY 4-Compartment Model", 11],
                       "simple_error.csv" => ["Harmonic Oscillator Model", 4],
                       "LV_error.csv" => ["Lotka-Volterra Model", 5])
model_legend_dict = Dict("FHN_error.csv" => :bottomleft,
                         "HIV_error.csv" => :bottomleft,
                         "crauste_error.csv" => :topright,
                         "daisy_ex3_error.csv" => :topright,
                         "daisy_mamil3_error.csv" => :bottomleft,
                         "daisy_mamil4_error.csv" => :topright,
                         "simple_error.csv" => :topright,
                         "LV_error.csv" => :topright)

for (amigo_f, pe_f, sciml_f, iqm_f) in zip(amigo, pe, sciml, iqm)
    println(amigo_f, " ", pe_f, " ", sciml_f)
    amigo_data = read_csv_data(amigo_dir * amigo_f)
    param_estim = read_txt_data(param_estim_dir * pe_f)
    sciml_data = read_txt_data(sciml_dir * sciml_f)
    iqm_data = read_csv_data(iqm_dir * iqm_f)
    title, num_params = model_name_dict[amigo_f]
    full_title = title * ", " * string(num_params) * " unknowns"
    scatter(amigo_data[!, "num_points"], amigo_data[!, "max_abs_rel_error"],
            xlabel = "Number of data points", ylabel = "Max. rel. err. [%] (log scale)",
            title = full_title, yaxis = :log, label = "AMIGO2", dpi = 300,
            markershape = :circle, color = :blue)
    scatter!(param_estim[!, "num_points"][1:19], param_estim[!, "max_abs_rel_err"][1:19],
             xlabel = "Number of data points", ylabel = "Max. rel. err. [%] (log scale)",
             title = full_title, label = "ParameterEstimation.jl",
             yaxis = :log, dpi = 300, markershape = :dtriangle, color = :red)
    scatter!(sciml_data[!, "num_points"], sciml_data[!, "max_abs_rel_err"],
             xlabel = "Number of data points", ylabel = "Max. rel. err. [%] (log scale)",
             title = full_title, label = "SciMLSensitivity.jl", yaxis = :log, dpi = 300,
             markershape = :square, color = :green)
    scatter!(iqm_data[!, "num_points"], iqm_data[!, "max_abs_rel_error"],
             xlabel = "Number of data points", ylabel = "Max. rel. err. [%] (log scale)",
             title = full_title, label = "IQM", yaxis = :log, dpi = 300,
             markershape = :diamond, color = :purple)
    plot!(amigo_data[!, "num_points"], amigo_data[!, "max_abs_rel_error"], yaxis = :log,
          label = "", dpi = 300, color = :blue)
    plot!(param_estim[!, "num_points"][1:19], param_estim[!, "max_abs_rel_err"][1:19],
          yaxis = :log, label = "", dpi = 300, color = :red)
    plot!(sciml_data[!, "num_points"], sciml_data[!, "max_abs_rel_err"], yaxis = :log,
          label = "", legend = model_legend_dict[amigo_f], dpi = 300, color = :green)
    plot!(iqm_data[!, "num_points"], iqm_data[!, "max_abs_rel_error"], yaxis = :log,
          label = "", dpi = 300, legend = model_legend_dict[amigo_f], color = :purple)
    # set y axis limits
    # ylims!(1e-4, 1e8)
    savefig(full_title * ".png")
end

# plot!(size = (1024, 768))