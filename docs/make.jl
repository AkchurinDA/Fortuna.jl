using Documenter
using DocumenterCitations
using Fortuna
using Distributions

Bibliography = CitationBibliography(
    joinpath(@__DIR__, "src/References.bib")
)

makedocs(
    sitename="Fortuna.jl",
    authors="Damir Akchurin, AkchurinDA@gmail.com",
    format=Documenter.HTML(
        sidebar_sitename=false,
        assets=String["assets/Citations.css"]
    ),
    pages=[
        "Home" => "index.md",
        "Random Variables" => "RandomVariables.md",
        "Isoprobabilistic Transformations" => [
            "Nataf Transformation" => "NatafTransformation.md",
            "Rosenblatt Transformation" => "RosenblattTransformation.md",
        ],
        "Reliability Analysis" => [
            "First-Order Reliability Methods (FORM)" => "FORM.md",
            "Second-Order Reliability Methods (SORM)" => "SORM.md",
        ],
        "Sensitivity Analysis" => "SensitivityAnalysis.md",
        "Examples" => "Examples.md",
        "References" => "References.md"
    ],
    plugins=[Bibliography]
)

