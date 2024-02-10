using Fortuna
using Documenter, DocumenterCitations

Bibliography = CitationBibliography(
    joinpath(@__DIR__, "src/References.bib"),
    style=:authoryear)

makedocs(
    sitename="Fortuna.jl",
    authors="Damir Akchurin, AkchurinDA@gmail.com",
    format=Documenter.HTML(
        sidebar_sitename=false,
        assets=[
            "assets/Favicon.ico",
            "assets/Custom.css"]),
    pages=[
        "Home" => "index.md",
        "Random Variables" => [
            "Generating Random Variables" => "RandomVariables/GenerateRV.md",
            "Sampling Random Variables" => "RandomVariables/SampleRV.md"
        ],
        "Isoprobabilistic Transformations" => [
            "Nataf Transformation" => "IsoprobabilisticTransformations/NatafTransformation.md",
            "Rosenblatt Transformation" => "IsoprobabilisticTransformations/RosenblattTransformation.md"
        ],
        "Reliability Problems" => [
            "Defining Reliability Problems" => "SolvingReliabilityProblems/DefiningReliabilityProblems.md",
            "Monte Carlo Simulations" => "SolvingReliabilityProblems/MCS.md",
            "First-Order Reliability Methods" => "SolvingReliabilityProblems/FORM.md",
            "Second-Order Reliability Methods" => "SolvingReliabilityProblems/SORM.md",
            "Subset Simulation Method" => "SolvingReliabilityProblems/SSM.md"
        ],
        "Inverse Reliability Problems" => [
            "Defining Inverse Reliability Problems" => "SolvingInverseReliabilityProblems/DefiningInverseReliabilityProblems.md"
        ],
        "Sensitivity Analysis" => "SensitivityAnalysis.md",
        "Examples" => "Examples.md",
        "Showcases" => "Showcases.md",
        "References" => "References.md"],
    plugins=[Bibliography])

deploydocs(
    repo="github.com/AkchurinDA/Fortuna.jl.git")