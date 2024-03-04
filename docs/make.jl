using Fortuna
using Documenter, DocumenterCitations

Bibliography = CitationBibliography(
    joinpath(@__DIR__, "src/References.bib"),
    style = :authoryear)

makedocs(
    sitename    = "Fortuna.jl",
    authors     = "Damir Akchurin, AkchurinDA@gmail.com",
    format      = Documenter.HTML(
        sidebar_sitename    = false,
        assets              = [
            "assets/Favicon.ico",
            "assets/Custom.css"
            ]),
    pages = [
        "Home"                                  => "index.md",
        "Random Variables"                      => [
            "Generating Random Variables"       => "RandomVariables/GenerateRV.md",
            "Sampling Random Variables"         => "RandomVariables/SampleRV.md"
            ],
        "Isoprobabilistic Transformations"      => [
            "Nataf Transformation"              => "IsoprobabilisticTransformations/NatafTransformation.md",
            "Rosenblatt Transformation"         => "IsoprobabilisticTransformations/RosenblattTransformation.md"
            ],
        "Reliability Problems"                  => [
            "Defining Reliability Problems"     => "ReliabilityProblems/DefiningReliabilityProblems.md",
            "Monte Carlo Simulations"           => "ReliabilityProblems/MCS.md",
            "Importance Sampling"               => "ReliabilityProblems/IS.md",
            "First-Order Reliability Methods"   => "ReliabilityProblems/FORM.md",
            "Second-Order Reliability Methods"  => "ReliabilityProblems/SORM.md",
            "Subset Simulation Method"          => "ReliabilityProblems/SSM.md"
            ],
        "Inverse Reliability Problems"          => "InverseReliabilityProblems.md",
        "Sensitivity Analysis"                  => "SensitivityAnalysis.md",
        "Examples"                              => "Examples.md",
        "Showcases"                             => "Showcases.md",
        "References"                            => "References.md"],
    plugins = [Bibliography])

deploydocs(
    repo = "github.com/AkchurinDA/Fortuna.jl.git")