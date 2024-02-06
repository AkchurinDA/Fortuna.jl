using Fortuna
using Documenter, DocumenterCitations

Bibliography = CitationBibliography(
    joinpath(@__DIR__, "src/References.bib"),
    style=:authoryear
)

makedocs(
    sitename="Fortuna.jl",
    authors="Damir Akchurin, AkchurinDA@gmail.com",
    format=Documenter.HTML(
        sidebar_sitename=false,
        assets=["assets/favicon.ico"]
    ),
    pages=[
        "Home" => "index.md",
        "Random Variables" => [
            "Generating Random Variables" => "GenerateRV.md",
            "Sampling Random Variables" => "SampleRV.md"
        ],
        "Isoprobabilistic Transformations" => [
            "Nataf Transformation" => "NatafTransformation.md",
            "Rosenblatt Transformation" => "RosenblattTransformation.md",
        ],
        "Reliability Analysis" => [
            "Overview" => "Overview.md",
            "First-Order Reliability Methods" => "FORM.md",
            "Second-Order Reliability Methods" => "SORM.md",
            "Subset Simulation Method" => "SSM.md"
        ],
        "Sensitivity Analysis" => "SensitivityAnalysis.md",
        "Examples" => "Examples.md",
        "Showcases" => "Showcases.md",
        "References" => "References.md"
    ],
    plugins=[Bibliography]
)

deploydocs(
    repo="github.com/AkchurinDA/Fortuna.jl.git"
)