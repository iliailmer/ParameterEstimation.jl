using Documenter, ParameterEstimation

makedocs(sitename = "ParameterEstimation.jl",
    modules = [ParameterEstimation],
    checkdocs = :exports,
    format = Documenter.HTML(edit_link = "main"),
    pages = ["Home" => "index.md",
        "Tutorial" => "tutorials/estimate.md",
        "API reference" => "api.md"])

deploydocs(repo = "github.com/iliailmer/ParameterEstimation.jl.git", devbranch = "main")
