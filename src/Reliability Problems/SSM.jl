"""
    SSM <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using Subset Simulation Method (SSM).

$(TYPEDFIELDS)
"""
Base.@kwdef struct SSM <: AbstractReliabililyAnalysisMethod
    "Probability of failure for each subset ``P_{0}``"
    P₀::Real = 0.10
    "Number of samples generated within each subset ``N``"
    NumSamples::Integer = 10 ^ 6
    "Maximum number of subsets ``M``"
    MaxNumSubsets::Integer = 50
end

"""
    SSMCache

Type used to perform reliability analysis using Subset Simulation Method (SSM).

$(TYPEDFIELDS)
"""
struct SSMCache
    "Samples generated within each subset in ``X``-space"
    XSamplesSubset::Vector{Matrix{Float64}}
    "Samples generated within each subset in ``U``-space"
    USamplesSubset::Vector{Matrix{Float64}}
    "Threshold for each subset ``C_{i}``"
    CSubset::Vector{Float64}
    "Probability of failure for each subset ``P_{f_{i}}``"
    PoFSubset::Vector{Float64}
    "Probability of failure ``P_{f}``"
    PoF::Float64
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::SSM)

Function used to solve reliability problems using Subset Simulation Method (SSM).
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::SSM)
    # Extract analysis details:
    P₀            = AnalysisMethod.P₀
    NumSamples    = AnalysisMethod.NumSamples
    MaxNumSubsets = AnalysisMethod.MaxNumSubsets

    # Extract problem data:
    X  = Problem.X
    ρˣ = Problem.ρˣ
    g  = Problem.g

    # Perform Nataf Transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Compute number of dimensions: 
    NumDimensions = length(X)

    # Compute number of Markov chains within each subset:
    NumMarkovChains = floor(Integer, P₀ * NumSamples)

    # Compute number of samples to generate within each Markov chain:
    NumSamplesChain = floor(Integer, NumSamples / NumMarkovChains)

    # Preallocate:
    USamplesSubset = Vector{Matrix{Float64}}()
    XSamplesSubset = Vector{Matrix{Float64}}()
    CSubset        = Vector{Float64}(undef, MaxNumSubsets)
    PoFSubset      = Vector{Float64}(undef, MaxNumSubsets)

    # Loop through each subset:
    for i in 1:MaxNumSubsets
        if i == 1
            # Generate samples in the standard normal space:
            USamples = Distributions.randn(NumDimensions, NumSamples)
        else
            # Preallocate:
            USamples = zeros(NumDimensions, NumMarkovChains * NumSamplesChain)

            # Generate samples using the Modified Metropolis-Hastings algorithm:
            for j in 1:NumMarkovChains
                USamples[:, (NumSamplesChain * (j - 1) + 1):(NumSamplesChain * j)] = 
                ModifiedMetropolisHastings(USamplesSubset[i - 1][:, j], CSubset[i - 1], NumDimensions, NumSamplesChain, NatafObject, g)
            end
        end

        # Evaluate the limit state function at the generated samples:
        GSamples = G(g, NatafObject, USamples)

        # Sort the values of the limit state function:
        GSamplesSorted = sort(GSamples)

        # Compute the threshold:
        CSubset[i] = Distributions.quantile(GSamplesSorted, P₀)

        # Check for convergance:
        if CSubset[i] < 0
            # Redefine the threshold:
            CSubset[i] = 0

            # Retain samples below the threshold:
            Indices = findall(x -> x ≤ CSubset[i], GSamples)
            push!(USamplesSubset, USamples[:, Indices])
            push!(XSamplesSubset, transformsamples(NatafObject, USamples[:, Indices], "U2X"))

            # Compute the probability of failure:
            PoFSubset[i] = length(Indices) / size(GSamples)[1]

            # Clean up the result:
            CSubset   = CSubset[1:i]
            PoFSubset = PoFSubset[1:i]

            # Break out:
            break
        else
            # Check for convergance:
            if i == MaxNumSubsets
                error("SSM did not converge. Try increasing the maximum number of subsets (MaxNumSubsets) or number of samples to generate within with subset (NumSamples).")
            end

            # Retain samples below the threshold:
            Indices = findall(x -> x ≤ CSubset[i], GSamples)
            push!(USamplesSubset, USamples[:, Indices])
            push!(XSamplesSubset, transformsamples(NatafObject, USamples[:, Indices], "U2X"))

            # Compute the probability of failure:
            PoFSubset[i] = length(Indices) / size(GSamples)[1]        
        end
    end

    # Compute the final probability of failure:
    PoF = prod(PoFSubset)

    # Return the result:
    return SSMCache(XSamplesSubset, USamplesSubset, CSubset, PoFSubset, PoF)
end

function ModifiedMetropolisHastings(StartingPoint::Vector{Float64}, CurrentThreshold::Float64, NumDimensions::Integer, NumSamplesChain::Integer, NatafObject::NatafTransformation, g::Function)
    # Preallocate:
    ChainSamples        = zeros(NumDimensions, NumSamplesChain)
    ChainSamples[:, 1]  = StartingPoint

    # Define a standard multivariate normal PDF:
    M = zeros(NumDimensions)
    Σ = Matrix(1.0 * LinearAlgebra.I, NumDimensions, NumDimensions)
    ϕ = Distributions.MvNormal(M, Σ)

    # Pregenerate uniformly-distributed samples:
    U = Distributions.rand(NumSamplesChain)

    # Generate samples:
    for i in 1:(NumSamplesChain - 1)
        # Define a proposal density:
        M = ChainSamples[:, i]
        Σ = Matrix(1.0 * LinearAlgebra.I, NumDimensions, NumDimensions)
        q = Distributions.MvNormal(M, Σ)

        # Propose a new state:
        ProposedState = Distributions.rand(q, 1)
        ProposedState = vec(ProposedState)

        # Compute the indicator function:
        XProposedState = transformsamples(NatafObject, ProposedState, "U2X")
        GProposedState = g(XProposedState)
        if GProposedState ≤ CurrentThreshold
            IF = 1
        else
            IF = 0
        end

        # Accept or reject the proposed state:
        α = (pdf(ϕ, ProposedState) * IF) / pdf(ϕ, ChainSamples[:, i]) # Acceptance ratio
        if U[i] <= α # Accept
            ChainSamples[:, i + 1] = ProposedState
        else # Reject
            ChainSamples[:, i + 1] = ChainSamples[:, i]
        end
    end

    # Return the result:
    return ChainSamples
end

function G(g::Function, NatafObject::NatafTransformation, USamples::AbstractMatrix)
    # Transform the samples:
    XSamples = transformsamples(NatafObject, USamples, "U2X")

    # Clean up the transformed samples:
    XSamplesClean = eachcol(XSamples)
    XSamplesClean = Vector.(XSamplesClean)
    
    # Evaluate the limit state function at the transform samples:
    GSamples = g.(XSamplesClean)

    # Return the result:
    return GSamples
end