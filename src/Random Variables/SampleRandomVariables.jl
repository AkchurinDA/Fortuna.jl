# Inverse Transform Sampling:
"""
    struct ITS <: AbstractSamplingTechnique

Type used to perform the Inverse Transform Sampling.
"""
struct ITS <: AbstractSamplingTechnique end

# Latin Hypercube Sampling:
"""
    struct LHS <: AbstractSamplingTechnique

Type used to perform the Latin Hypercube Sampling.
"""
struct LHS <: AbstractSamplingTechnique end

# --------------------------------------------------
# GENERATE SAMPLES FROM A RANDOM VARIABLE
# --------------------------------------------------
"""
    rand(RNG::Distributions.AbstractRNG, RandomVariable::Distributions.UnivariateDistribution, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)

Function used to generate samples from an *random variable*.
If `SamplingTechnique` is:
- `ITS()` samples are generated using Inverse Transform Sampling technique.
- `LHS()` samples are generated using Latin Hypercube Sampling technique.
"""
function Distributions.rand(RNG::Distributions.AbstractRNG, RandomVariable::Distributions.UnivariateDistribution, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
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
        Samples = Distributions.quantile.(RandomVariable, UniformSamples)
    end

    # Return the result:
    return Samples
end

Distributions.rand(RandomVariable::Distributions. UnivariateDistribution, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique) = 
Distributions.rand(Distributions.default_rng(), RandomVariable, NumSamples, SamplingTechnique)

# --------------------------------------------------
# GENERATE SAMPLES FROM A RANDOM VECTOR
# --------------------------------------------------
"""
    rand(RNG::Distributions.AbstractRNG, RandomVector::Vector{<:Distributions.UnivariateDistribution}, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)

Function used to generate samples from a *random vector with uncorrelated marginals*.
If `SamplingTechnique` is:
- `ITS()` samples are generated using Inverse Transform Sampling technique.
- `LHS()` samples are generated using Latin Hypercube Sampling technique.
"""
function Distributions.rand(RNG::Distributions.AbstractRNG, RandomVector::Vector{<:Distributions.UnivariateDistribution}, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
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

Distributions.rand(RandomVector::Vector{<:Distributions. UnivariateDistribution}, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique) = 
Distributions.rand(Distributions.default_rng(), RandomVector, NumSamples, SamplingTechnique)

# --------------------------------------------------
# GENERATE SAMPLES FROM A TRANSFORMATION OBJECT
# --------------------------------------------------
# Nataf Transformation:
"""
    rand(RNG::Distributions.AbstractRNG, TransformationObject::NatafTransformation, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)

Function used to generate samples from a *random vector with correlated marginals* using Nataf Transformation object.
If `SamplingTechnique` is:
- `ITS()` samples are generated using Inverse Transform Sampling technique.
- `LHS()` samples are generated using Latin Hypercube Sampling technique.
"""
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

Distributions.rand(TransformationObject::NatafTransformation, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique) = 
Distributions.rand(Distributions.default_rng(), TransformationObject::NatafTransformation, NumSamples::Int, SamplingTechnique::AbstractSamplingTechnique)
