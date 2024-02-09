# Subset Simluation Method:
"""
    analyze(Problem::ReliabilityProblem, AnalysisMethod::SSM)

The function estimates small probabilities of failure using Subset Simulation Method (SSM).
"""
function analyze(Problem::ReliabilityProblem, AnalysisMethod::SSM)
    # Extract the analysis method:
    P₀ = AnalysisMethod.P₀
    NumSamples = AnalysisMethod.NumSamples
    MaxNumSubsets = AnalysisMethod.MaxNumSubsets

    # Extract the problem data:
    X = Problem.X
    ρˣ = Problem.ρˣ
    g = Problem.g

    # Compute the number of marginal distributions:
    NumDims = length(X)

    # Compute the number of samples to keep within each subset:
    NumSamplesKeep = floor(Integer, P₀ * NumSamples)

    # Compute the number of samples to generate within each chain:
    NumSamplesChain = floor(Integer, NumSamples / NumSamplesKeep)

    # Preallocate:
    USamplesSubset  = Vector{Matrix{Float64}}()
    XSamplesSubset  = Vector{Matrix{Float64}}()
    CSubset         = Vector{Float64}(undef, MaxNumSubsets)
    PoFSubset       = Vector{Float64}(undef, MaxNumSubsets)

    # Perform the Nataf Transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Loop through each subset:
    for i in 1:MaxNumSubsets
        if i == 1
            # Generate samples in the standard normal space:
            USamples = randn(NumSamples, NumDims)
        else
            # Preallocate:
            USamples = zeros(NumSamplesKeep * NumSamplesChain, NumDims)

            # Generate MCMCs according to the modified Metropolis-Hastings algorithm:
            for j in 1:NumSamplesKeep
                USamples[(NumSamplesChain*(j-1)+1):(NumSamplesChain*j), :] = MMH(USamplesSubset[i-1][j, :], CSubset[i-1], NumDims, NumSamplesChain, NatafObject, g)
            end
        end

        # Evaluate the limit state function at the generated samples:
        GSamples = G(g, NatafObject, USamples)

        # Sort the values of the limit state function:
        GSamplesSorted = sort(GSamples)

        # Compute the threshold:
        CSubset[i] = quantile(GSamplesSorted, P₀)

        # Check for convergance:
        if CSubset[i] > 0
            # Retain samples below the threshold:
            Indices = findall(x -> x ≤ CSubset[i], GSamples)
            push!(USamplesSubset, USamples[Indices, :])
            push!(XSamplesSubset, transformsamples(NatafObject, USamples[Indices, :], "U2X"))

            # Compute the probability of failure:
            PoFSubset[i] = length(Indices) / size(GSamples)[1]
        else
            # Redefine the threshold:
            CSubset[i] = 0

            # Retain samples below the threshold:
            Indices = findall(x -> x ≤ CSubset[i], GSamples)
            push!(USamplesSubset, USamples[Indices, :])
            push!(XSamplesSubset, transformsamples(NatafObject, USamples[Indices, :], "U2X"))

            # Compute the probability of failure:
            PoFSubset[i] = length(Indices) / size(GSamples)[1]

            # Clean up the result:
            CSubset = CSubset[1:i]
            PoFSubset = PoFSubset[1:i]

            # Break out:
            break
        end
    end

    # Compute the probability of failure:
    PoF = prod(PoFSubset)

    # Return the result:
    return SSMCache(XSamplesSubset, USamplesSubset, CSubset, PoFSubset, PoF)
end

function MMH(StartingPoint::Vector{Float64}, CurrentThreshold::Float64, NumDims::Integer, NumSamples::Integer, NatafObject::NatafTransformation, g::Function)
    # Preallocate:
    ChainSamples        = zeros(NumSamples, NumDims)
    ChainSamples[1, :]  = StartingPoint

    # Define a standard multivariate normal PDF:
    M       = zeros(NumDims)
    Σ       = I(NumDims)
    MVN     = MvNormal(M, Σ)

    # Pregenerate uniformly-distributed samples:
    U = rand(NumSamples)

    # Generate samples:
    for i in 1:NumSamples-1
        # Define a proposal density:
        M = ChainSamples[i, :]
        Σ = I(NumDims)
        ProposalDensity = MvNormal(M, Σ)

        # Propose a new state:
        ProposedState = rand(ProposalDensity, 1)
        ProposedState = vec(ProposedState)

        # Compute the indicator function:
        XProposedState = transformsamples(NatafObject, ProposedState, "U2X")
        gProposedState = g(XProposedState)
        if gProposedState ≤ CurrentThreshold
            IF = 1
        else
            IF = 0
        end

        # Accept or reject the proposed state:
        α = (pdf(MVN, ProposedState) * IF) / pdf(MVN, ChainSamples[i, :]) # Acceptance ratio
        if U[i] <= α # Accept
            ChainSamples[i+1, :] = ProposedState
        else # Reject
            ChainSamples[i+1, :] = ChainSamples[i, :]
        end
    end

    # Return the result:
    return ChainSamples
end

function G(g::Function, NatafObject::NatafTransformation, USamples::Matrix{Float64})
    # Transform samples:
    XSamples = transformsamples(NatafObject, USamples, "U2X")

    # Clean up the transformed samples:
    XSamples = eachrow(XSamples)
    XSamples = Vector.(XSamples)

    # Evaluate the limit state function at the transform samples:
    GSamples = g.(XSamples)

    # Return the result:
    return GSamples
end