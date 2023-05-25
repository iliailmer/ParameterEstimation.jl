import json


class Benchmarks:
    def __init__(self, models_filename):
        self.models = []
        with open(models_filename) as models_json_file:
            models_json_obj = json.load(models_json_file)
            for model_obj in models_json_obj["models"]:
                self.models.append(Model(**model_obj))


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
        pass


if __name__ == "__main__":
    b = Benchmarks("models.json")
