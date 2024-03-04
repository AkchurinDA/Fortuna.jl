"""
    analyze(Problem::ReliabilityProblem, AnalysisMethod::IS)

The function solves the provided reliability problem using Importance Sampling technique.
"""
function analyze(Problem::ReliabilityProblem, AnalysisMethod::IS)
    # Extract the analysis details:
    q                   = AnalysisMethod.q
    NumSamples          = AnalysisMethod.NumSamples

    # Extract data:
    g   = Problem.g
    X   = Problem.X
    ρˣ  = Problem.ρˣ

    # Check dimensions:
    if length(q) != length(X)
        error("Dimensionality of the proposal distribution does not match the dimensionality of the random vector.")
    end

    # If the marginal distrbutions are correlated, define a Nataf object:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    XSamples = rand(q, NumSamples)

    # Clean up the generated samples:
    XSamples        = transpose(XSamples)
    XSamples        = Matrix(XSamples)
    XSamplesClean   = eachrow(XSamples)
    XSamplesClean   = Vector.(XSamplesClean)

    # Evaluate the target and proposal probability density functions at the generate samples:
    fSamples = jointpdf.(NatafObject, XSamplesClean)
    qSamples = pdf(q, XSamplesClean)

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