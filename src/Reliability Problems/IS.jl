"""
    IS <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using Importance Sampling (IS) method.
"""
Base.@kwdef struct IS <: AbstractReliabililyAnalysisMethod
    "Proposal probability density function ``q``"
    q::Distributions.ContinuousMultivariateDistribution
    "Number of samples to generate ``N``"
    NumSimulations::Integer = 1E6
end

"""
    ISCache

Type used to store results of reliability analysis performed using Importance Sampling (IS) method.
"""
struct ISCache
    "Samples generated from the proposal probability density function ``\\vec{x}_{i}``"
    Samples::Matrix{Float64}
    "Target probability density function evaluated at each sample ``f(\\vec{x}_{i})``"
    fValues::Vector{Float64}
    "Proposal probability density function evaluated at each sample ``q(\\vec{x}_{i})``"
    qValues::Vector{Float64}
    "Limit state function evalued at each sample ``g(\\vec{x}_{i})``"
    gValues::Vector{Float64}
    "Probability of failure ``P_{f}``"
    PoF::Float64
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::IS; showprogressbar = false)

Function used to solve reliability problems using Importance Sampling (IS) method.
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::IS; showprogressbar = false)
    # Extract the analysis details:
    q              = AnalysisMethod.q
    NumSimulations = AnalysisMethod.NumSimulations

    # Extract data:
    g  = Problem.g
    X  = Problem.X
    ρˣ = Problem.ρˣ

    # Error-catching:
    length(q) == length(X) || throw(DimensionMismatch("Dimensionality of the proposal distribution does not match the dimensionality of the random vector!"))

    # If the marginal distrbutions are correlated, define a Nataf object:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    Samples = rand(q, NumSimulations)

    # Evaluate the target and proposal probability density functions at the generate samples:
    fValues = pdf(NatafObject, Samples)
    qValues = pdf(q, Samples)

    # Evaluate the limit state function at the generate samples:
    gValues = Vector{Float64}(undef, NumSimulations)
    ProgressMeter.@showprogress desc = "Evaluating the limit state function..." enabled = showprogressbar for i in axes(Samples, 2)
        gValues[i] = g(Samples[:, i])
    end

    # Evaluate the indicator function at the generate samples:
    IValues = gValues .≤ 0
    IValues = Int.(IValues)

    # Compute the probability of failure:
    PoF = mean((IValues .* fValues) ./ qValues)

    # Return results:
    return ISCache(Samples, fValues, qValues, gValues, PoF)
end