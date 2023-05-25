import json
import csv
import os
data_dir = "data"

class Benchmarks:
    def __init__(self, models_filename):
        self.models = []
        with open(models_filename) as models_json_file:
            models_json_obj = json.load(models_json_file)
            for model_obj in models_json_obj["models"]:
                self.models.append(Model(**model_obj))

    def run_amigo(self):
        from amigo import run
        for model in self.models:
            run(model)

    def run_iqm(self):
        from iqm import run
        for model in self.models:
            run(model)

    def run(self):
        self.run_amigo()
        self.run_iqm()

class Model:
    def __init__(
        self,
        name,
        start,
        end,
        interval_count,
        state_variables,
        parameter_variables,
        diff_eqns,
        output_variables,
        output_eqns,
        output_data_filename,
    ):
        self.output_data_filename = output_data_filename
        self.__read_data_file()

    def __read_data_file(self):
        with open(os.path.join(data_dir, self.output_data_filename), newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                assert 't' in row
                for varname in self.output_variables:
                    assert f'{varname}(t)' in row
                print(row)


if __name__ == "__main__":
    b = Benchmarks("models.json")
    b.run()
