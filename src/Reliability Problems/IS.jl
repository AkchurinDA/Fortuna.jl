"""
    IS <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using Importance Sampling (IS) method.
"""
Base.@kwdef struct IS <: AbstractReliabililyAnalysisMethod
    "Proposal probability density function ``q``"
    q::Distributions.Sampleable
    "Number of samples to generate ``N``"
    NumSamples::Integer = 10 ^ 6
end

"""
    ISCache

Type used to store results of reliability analysis performed using Importance Sampling (IS) method.
"""
struct ISCache
    "Generated samples"
    Samples::Matrix{Float64}
    "Probability of failure ``P_{f}``"
    PoF::Float64
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::IS)

Function used to solve reliability problems using Importance Sampling (IS) method.
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::IS)
    # Extract the analysis details:
    q          = AnalysisMethod.q
    NumSamples = AnalysisMethod.NumSamples

    # Extract data:
    g  = Problem.g
    X  = Problem.X
    ρˣ = Problem.ρˣ

    # Check dimensions:
    if length(q) != length(X)
        error("Dimensionality of the proposal distribution does not match the dimensionality of the random vector.")
    end

    # If the marginal distrbutions are correlated, define a Nataf object:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    XSamples = rand(q, NumSamples)

    # Clean up the generated samples:
    XSamplesClean = eachcol(XSamples)
    XSamplesClean = Vector.(XSamplesClean)

    # Evaluate the target and proposal probability density functions at the generate samples:
    fSamples = pdf(NatafObject, XSamples)
    qSamples = pdf(q, XSamples)

    # Evaluate the limit state function at the generate samples:
    gSamples = g.(XSamplesClean)

    # Evaluate the indicator function at the generate samples:
    ISamples = gSamples .≤ 0
    ISamples = Int.(ISamples)

    # Compute the probability of failure:
    PoF = mean((ISamples .* fSamples) ./ qSamples)

    # Return results:
    return ISCache(XSamples, PoF)
end