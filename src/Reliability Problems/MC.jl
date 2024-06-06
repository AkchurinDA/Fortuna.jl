"""
    MC <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using Monte Carlo (MC) simulations.

$(TYPEDFIELDS)
"""
Base.@kwdef struct MC <: AbstractReliabililyAnalysisMethod
    "Number of simulations ``N``"
    NumSimulations::Integer = 1E6
    "Sampling technique"
    SamplingTechnique::Symbol = :LHS
end

"""
    MCCache

Type used to store results of reliability analysis performed using Monte Carlo (MC) simulations.

$(TYPEDFIELDS)
"""
struct MCCache
    "Generated samples ``\\vec{x}_{i}``"
    Samples::Matrix{Float64}
    "Limit state function evalued at each sample ``g(\\vec{x}_{i})``"
    gValues::Vector{Float64}
    "Probability of failure ``P_{f}``"
    PoF::Float64
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::MC; showprogressbar = false)

Function used to solve reliability problems using Monte Carlo (MC) simulations.
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::MC; showprogressbar = false)
    # Extract the analysis details:
    NumSimulations    = AnalysisMethod.NumSimulations
    SamplingTechnique = AnalysisMethod.SamplingTechnique

    # Extract data:
    g  = Problem.g
    X  = Problem.X
    ρˣ = Problem.ρˣ

    # If the marginal distrbutions are correlated, define a Nataf object:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    if (SamplingTechnique != :ITS) && (SamplingTechnique != :LHS)
        throw(ArgumentError("Provided sampling technique is not supported!"))
    else
        Samples, _, _ = rand(NatafObject, NumSimulations, SamplingTechnique)
    end

    # Evaluate the limit state function at the generate samples:
    gValues = Vector{Float64}(undef, NumSimulations)
    ProgressMeter.@showprogress desc = "Running Monte Carlo simulations..." enabled = showprogressbar for i in axes(Samples, 2)
        gValues[i] = g(Samples[:, i])
    end

    # Compute the probability of failure:
    PoF = count(x -> x ≤ 0, gValues) / NumSimulations

    # Return results:
    return MCCache(Samples, gValues, PoF)
end