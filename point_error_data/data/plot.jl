amigo_dir = "./amigo/"
param_estim_dir = "./pe.jl/"
sciml_dir = "./sciml/"
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

for (amigo_f, pe_f, sciml_f) in zip(amigo, pe, sciml)
    amigo_data = read_csv_data(amigo_dir * amigo_f)
    param_estim = read_txt_data(param_estim_dir * pe_f)
    sciml_data = read_txt_data(sciml_dir * sciml_f)
    break
end