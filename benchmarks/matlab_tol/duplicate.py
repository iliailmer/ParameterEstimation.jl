import os
from glob import glob

import pandas as pd


def dup_dir(src, max_idx):
    idx = 5
    while idx < max_idx:
        dst = src.split("_")[0] + "_" + str(idx)
        os.system("cp -r " + src + " " + dst)
        idx += 1


def convert_to_csv(dir):
    files = glob(os.path.join(dir, "*.txt"))
    for file in files:
        with open(file, "r") as f:
            lines = f.readlines()
            lines = [line.split() for line in lines]
            lines = [",".join(line) for line in lines]
            with open(file[:-4] + ".csv", "w") as f:
                f.write("\n".join(lines))


def rearrange_cols(file):
    df = pd.read_csv(file)
    df = df[["t"] + [col for col in df.columns if col != "t"]]
    df.to_csv(file, index=False)


def edit_experiment_file(file, df):
    with open(file, "r") as f:
        data = f.read()
        # replace everything after [Values] with content of df
        tmp = df.to_string(index=False, header=False, float_format="%.16f").split("\n")
        tmp = "\n".join([",".join(x.split()) for x in tmp])
        data = data.split("[Values]")[0] + "[Values]\n" + tmp
        with open(file, "w") as f:
            f.write(data)


dirs = os.listdir(
    "/Users/iliailmer/projects/julia/ParameterEstimation.jl/sharable/matlab/point_error_data/samples"
)
csvs = glob(
    "/Users/iliailmer/projects/julia/ParameterEstimation.jl/sharable/matlab/point_error_data/samples/*/*.csv"
)
experiments = glob(
    "/Users/iliailmer/projects/julia/ParameterEstimation.jl/sharable/matlab/iqm/*/project/experiments/*/*.csv"
)

path_to_iqm_models = sorted(
    glob("/Users/iliailmer/projects/julia/ParameterEstimation.jl/sharable/matlab/iqm/*")
)
path_to_data_samples = sorted(
    glob(
        "/Users/iliailmer/projects/julia/ParameterEstimation.jl/sharable/matlab/point_error_data/samples/*"
    )
)
assert all(
    [
        (p1.split("/")[-1] == p2.split("/")[-1])
        for p1, p2 in zip(path_to_iqm_models, path_to_data_samples)
    ]
)
for model, sample in zip(path_to_iqm_models, path_to_data_samples):
    csv_samples = sorted(
        glob(os.path.join(sample, "*.csv")),
        key=lambda x: int(x.split("-")[-1].split(".")[0]),
    )
    experiments = sorted(
        glob(os.path.join(model, "project/experiments/*/*.csv")),
        key=lambda x: int(x.split("/")[-2].split("_")[-1]),
    )
    for experiment, sample in zip(experiments, csv_samples):
        df = pd.read_csv(sample)
        edit_experiment_file(experiment, df)
        print(experiment, sample)
        # break
    # break
