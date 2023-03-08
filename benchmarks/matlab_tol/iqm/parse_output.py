from collections import defaultdict
from glob import glob

import pandas as pd
import regex as re
from tqdm.auto import tqdm


def parse_output(output):
    """
    Match all lines that look like this
        Ab_1 : 1.2345e-01
    """
    res = re.findall(r"([A-Za-z0-9_]+)\s*=\s*([0-9.\+\-e]+)", output)
    #
    return res


output_files = glob("outputs/*.out")
pbar = tqdm(total=len(output_files))
for out in output_files:
    pbar.update(1)
    pbar.set_description(out)

    df = pd.DataFrame()  # create a dataframe to store the results
    for idx, num_points in enumerate(range(3, 22)):  # loop over the number of points
        df.loc[
            idx, "num_points"
        ] = num_points  # store the number of points in the dataframe
    results_dict = defaultdict(list)
    with open(out, "r") as f:
        output = f.read()
        res = parse_output(output)
        for each in res:
            results_dict[each[0]].append(float(each[1]))
    for key, value in results_dict.items():
        df[key] = value

    df.set_index("num_points", inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    df.to_csv(out.replace(".out", ".csv"), index=True)
    # break

# print(df)


for each in output_files:
    csv_file = pd.read_csv(each.replace(".out", ".csv"))
    truth_file = pd.read_csv(each.replace(".out", "_true.csv"))
    df = pd.DataFrame()  # create a dataframe to store the results
    df["num_points"] = range(3, 22)  # store the number of points in the dataframe
    truth = truth_file.iloc[0, 1:]
    for idx, row in csv_file.iterrows():
        df.loc[idx, "num_points"] = row["num_points"]
        row = row[1:]
        df.loc[idx, "max_abs_rel_error"] = max(
            [
                abs((truth[c] - row[c]) / truth[c]) * 100 for c in row.index
            ]  # error!!!!! check colunns
        )

    df.to_csv(each.replace(".out", "_error.csv"), index=False)
