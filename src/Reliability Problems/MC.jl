"""
    struct MC <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using Monte Carlo simulations.

$(TYPEDFIELDS)
"""
Base.@kwdef struct MC <: AbstractReliabililyAnalysisMethod
    "Number of samples ``N``"
    NumSamples          ::Integer = 10 ^ 6
    "Sampling technique"
    SamplingTechnique   ::AbstractSamplingTechnique = ITS()
end

"""
    struct MCCache

Type used to store results of reliability analysis performed using Monte Carlo simulations.

$(TYPEDFIELDS)
"""
struct MCCache
    "Generated samples"
    Samples ::Matrix{Float64}
    "Probability of failure ``P_{f}``"
    PoF     ::Float64
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::MC)

Function used to solve reliability analysis using Monte Carlo simulations.
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::MC)
    # Extract the analysis details:
    NumSamples          = AnalysisMethod.NumSamples
    SamplingTechnique   = AnalysisMethod.SamplingTechnique

    # Extract data:
    g   = Problem.g
    X   = Problem.X
    ρˣ  = Problem.ρˣ

    # If the marginal distrbutions are correlated, define a Nataf object:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    if !isa(SamplingTechnique, ITS) && !isa(SamplingTechnique, LHS)
        error("Invalid sampling technique.")
    else
        XSamples, _, _ = rand(NatafObject, NumSamples, SamplingTechnique)
    end

    # Clean up the generated samples:
    XSamplesClean = eachcol(XSamples)
    XSamplesClean = Vector.(XSamplesClean)

    # Evaluate the limit state function at the generate samples:
    gSamples = g.(XSamplesClean)

    # Compute the probability of failure:
    PoF = count(x -> x ≤ 0, gSamples) / NumSamples

    # Return results:
    return MCCache(XSamples, PoF)
end