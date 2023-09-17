using Documenter
using Fortuna
using Distributions

makedocs(
    sitename="Fortuna.jl",
    pages=[
        "Home" => "index.md",
        "Nataf Transformation" => "nataf.md",
        "Rosenblatt Transformation" => "rosenblatt.md"]
)