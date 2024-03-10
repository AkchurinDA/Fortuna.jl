function Distributions.rand(RNG::Distributions.AbstractRNG, RandomVariable::Distributions.Sampleable, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
    # Sample:
    if !isa(SamplingTechnique, ITS) && !isa(SamplingTechnique, LHS)
        error("Provided sampling technique is not supported.")
    elseif isa(SamplingTechnique, ITS)
        # Generate samples:
        Samples = Distributions.rand(RNG, RandomVariable, NumSamples)
    elseif isa(SamplingTechnique, LHS)
        # Define the lower limits of each strata:
        LowerLimits = collect(0:(1 / NumSamples):((NumSamples - 1) / NumSamples))

        # Generate samples within each strata:
        UniformSamples = LowerLimits + Distributions.rand(RNG, Distributions.Uniform(0, 1 / NumSamples), NumSamples)

        # Shuffle samples:
        UniformSamples = Random.shuffle(RNG, UniformSamples)

        # Generate samples:
        Samples = Distributions.quantile(RandomVariable, UniformSamples)
    end

    # Return the result:
    return Samples
end

Distributions.rand(RandomVariable::Distributions.Sampleable, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique) = Distributions.rand(Distributions.default_rng(), RandomVariable, NumSamples, SamplingTechnique)

# Generate samples from a random vector:
function Distributions.rand(RNG::Distributions.AbstractRNG, RandomVector::Vector{<:Distributions.Sampleable}, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
    # Compute number of dimensions:
    NumDimensions = length(RandomVector)

    # Preallocate:
    Samples = Matrix{Float64}(undef, NumDimensions, NumSamples)

    # Sample:
    for i in 1:NumDimensions
        Samples[i, :] = Distributions.rand(RNG, RandomVector[i], NumSamples, SamplingTechnique)
    end

    # Return the result:
    return Samples
end

Distributions.rand(RandomVector::Vector{<:Distributions.Sampleable}, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique) = Distributions.rand(Distributions.default_rng(), RandomVector, NumSamples, SamplingTechnique)

# Generate samples from a transformation object:
function Distributions.rand(RNG::Distributions.AbstractRNG, TransformationObject::NatafTransformation, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
    # Extract data:
    X = TransformationObject.X
    L = TransformationObject.L

    # Compute number of dimensions:
    NumDimensions = length(X)

    # Generate samples of uncorrelated normal random variables U:
    USamples = Matrix{Float64}(undef, NumDimensions, NumSamples)
    for i in 1:NumDimensions
        USamples[i, :] = Distributions.rand(RNG, Distributions.Normal(), NumSamples, SamplingTechnique)
    end

    # Generate samples of correlated normal random variables Z:
    ZSamples = L * USamples

    # Generate samples of correlated non-normal random variables X:
    XSamples = Matrix{Float64}(undef, NumDimensions, NumSamples)
    for i in 1:NumDimensions
        XSamples[i, :] = Distributions.quantile.(X[i], Distributions.cdf.(Distributions.Normal(), ZSamples[i, :]))
    end

    return XSamples, ZSamples, USamples
end

Distributions.rand(TransformationObject::NatafTransformation, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique) = Distributions.rand(Distributions.default_rng(), TransformationObject::NatafTransformation, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
