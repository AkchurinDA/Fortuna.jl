using Fortuna
using Documenter, DocumenterCitations
import Distributions

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
            "assets/Custom.css"
        ]),
    pages=[
        "Home"                                      => "index.md",
        "Random Variables"                          => [
            "Generating Random Variables"           => "Random Variables/GenerateRandomVariables.md",
            "Sampling Random Variables"             => "Random Variables/SampleRandomVariables.md"
        ],
        "Isoprobabilistic Transformations"          => [
            "Nataf Transformation"                  => "Isoprobabilistic Transformations/NatafTransformation.md",
            "Rosenblatt Transformation"             => "Isoprobabilistic Transformations/RosenblattTransformation.md"
        ],
        "Reliability Problems"                      => [
            "Reliability Problems"                  => "Reliability Problems/ReliabilityProblems.md",
            "Monte Carlo Methods"                   => [
                "Direct Monte Carlo Simulations"    => "Reliability Problems/MC.md",
                "Importance Sampling"               => "Reliability Problems/IS.md"
            ],
            "First-Order Reliability Methods"       => "Reliability Problems/FORM.md",
            "Second-Order Reliability Methods"      => "Reliability Problems/SORM.md",
            "Subset Simulation Method"              => "Reliability Problems/SSM.md"
        ],
        "Inverse Reliability Problems"              => "Reliability Problems/InverseReliabilityProblems.md",
        "Sensitivity Problems"                      => "Reliability Problems/SensitivityProblems.md",
        "Examples"                                  => "Examples.md",
        "Showcases"                                 => "Showcases.md",
        "References"                                => "References.md"
    ],
    plugins=[Bibliography])

deploydocs(
    repo="github.com/AkchurinDA/Fortuna.jl.git")