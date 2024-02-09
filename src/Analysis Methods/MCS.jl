# First-Order Reliability Method:
"""
    analyze(Problem::ReliabilityProblem, AnalysisMethod::MC)

The function solves the provided reliability problem using Monte Carlo simulations.
"""
function analyze(Problem::ReliabilityProblem, AnalysisMethod::MCS)
    # Extract the analysis details:
    NumSamples = AnalysisMethod.NumSamples
    SamplingTechnique = AnalysisMethod.SamplingTechnique

    # Extract data:
    g = Problem.g
    X = Problem.X
    ρˣ = Problem.ρˣ

    # If the marginal distrbutions are correlated, define a Nataf object:
    CorrelationFlag = false
    if any((ρˣ - I) .!= 0)
        NatafObject = NatafTransformation(X, ρˣ)
        CorrelationFlag = true
    end

    # Generate samples:
    if !isa(SamplingTechnique, ITS) && !isa(SamplingTechnique, LHS)
        error("Invalid sampling technique.")
    else
        if CorrelationFlag
            XSamples, _, _ = samplerv(NatafObject, NumSamples, SamplingTechnique)
        else
            XSamples = samplerv(X, NumSamples, SamplingTechnique)
        end
    end

    # Clean up the generated samples:
    XSamples = eachrow(XSamples)
    XSamples = Vector.(XSamples)

    # Evaluate the limit state function at the generate samples:
    gSamples = g.(XSamples)

    # Compute the probability of failure:
    PoF = count(x -> x ≤ 0, gSamples) / NumSamples

    # Return results:
    return MCSCache(PoF)
end