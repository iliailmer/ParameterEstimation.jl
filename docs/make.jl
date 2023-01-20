using Documenter, ParameterEstimation

makedocs(sitename = "ParameterEstimation.jl",
         pages = ["Home" => "index.md",
             "Tutorials" => [
                 "Estimation" => [
                     "Estimation for a single degree" => "tutorials/estimate_single_degree.md",
                     "Estimation for a range of degrees" => "tutorials/estimate.md",
                 ],
                 "Identifiability" => "tutorials/identifiability_tutorial.md",
             ],
             "Library" => [
                 "Estimation" => "library/estimate.md",
                 "Filtering" => "library/filtering.md",
                 "Identifiability" => "library/identifiability/identifiability.md",
                 "Rational Interpolation" => [
                     "Rational Interpolation" => "library/rational_interpolation/rational_interpolation.md",
                     "Interpolant" => "library/rational_interpolation/interpolant.md",
                 ],
             ]])