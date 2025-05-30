# # Inverse Transform Sampling:
# """
#     struct ITS <: AbstractSamplingTechnique

# Type used to perform the Inverse Transform Sampling.
# """
# struct ITS <: AbstractSamplingTechnique end

# # Latin Hypercube Sampling:
# """
#     struct LHS <: AbstractSamplingTechnique

# Type used to perform the Latin Hypercube Sampling.
# """
# struct LHS <: AbstractSamplingTechnique end

# --------------------------------------------------
# GENERATE SAMPLES FROM A RANDOM VARIABLE
# --------------------------------------------------
"""
    rand(RNG::Distributions.AbstractRNG, RandomVariable::Distributions.UnivariateDistribution, NumSamples::Integer, SamplingTechnique::Symbol)

Function used to generate samples from an *random variable*.
If `SamplingTechnique` is:
- `:ITS` samples are generated using Inverse Transform Sampling technique.
- `:LHS` samples are generated using Latin Hypercube Sampling technique.
"""
function Distributions.rand(RNG::Distributions.AbstractRNG, RandomVariable::Distributions.ContinuousUnivariateDistribution, NumSamples::Integer, SamplingTechnique::Symbol)
    # Generate samples:
    if (SamplingTechnique != :ITS) && (SamplingTechnique != :LHS)
        throw(ArgumentError("Provided sampling technique is not supported!"))
    elseif SamplingTechnique == :ITS
        # Generate samples:
        Samples = Distributions.rand(RNG, RandomVariable, NumSamples)
    elseif SamplingTechnique == :LHS
        # Define the lower limits of each strata:
        LowerLimits = collect(range(0, (NumSamples - 1) / NumSamples, NumSamples))

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

Distributions.rand(RandomVariable::Distributions.ContinuousUnivariateDistribution, NumSamples::Integer, SamplingTechnique::Symbol) = 
Distributions.rand(Distributions.default_rng(), RandomVariable, NumSamples, SamplingTechnique)

# --------------------------------------------------
# GENERATE SAMPLES FROM A RANDOM VECTOR
# --------------------------------------------------
"""
    rand(RNG::Distributions.AbstractRNG, RandomVector::Vector{<:Distributions.UnivariateDistribution}, NumSamples::Integer, SamplingTechnique::Symbol)

Function used to generate samples from a *random vector with uncorrelated marginals*.
If `SamplingTechnique` is:
- `:ITS` samples are generated using Inverse Transform Sampling technique.
- `:LHS` samples are generated using Latin Hypercube Sampling technique.
"""
function Distributions.rand(RNG::Distributions.AbstractRNG, RandomVector::Vector{<:Distributions.ContinuousUnivariateDistribution}, NumSamples::Integer, SamplingTechnique::Symbol)
    # Compute number of dimensions:
    NumDimensions = length(RandomVector)

    # Generate samples:
    Samples = [Distributions.rand(RNG, RandomVector[i], NumSamples, SamplingTechnique) for i in 1:NumDimensions]
    Samples = vcat(Samples'...)

    # Return the result:
    return Samples
end

Distributions.rand(RandomVector::Vector{<:Distributions.ContinuousUnivariateDistribution}, NumSamples::Integer, SamplingTechnique::Symbol) = 
Distributions.rand(Distributions.default_rng(), RandomVector, NumSamples, SamplingTechnique)

# --------------------------------------------------
# GENERATE SAMPLES FROM A TRANSFORMATION OBJECT
# --------------------------------------------------
# Nataf Transformation:
"""
    rand(RNG::Distributions.AbstractRNG, TransformationObject::NatafTransformation, NumSamples::Integer, SamplingTechnique::Symbol)

Function used to generate samples from a *random vector with correlated marginals* using Nataf Transformation object.
If `SamplingTechnique` is:
- `:ITS` samples are generated using Inverse Transform Sampling technique.
- `:LHS` samples are generated using Latin Hypercube Sampling technique.
"""
function Distributions.rand(RNG::Distributions.AbstractRNG, TransformationObject::NatafTransformation, NumSamples::Integer, SamplingTechnique::Symbol)
    # Extract data:
    X = TransformationObject.X
    L = TransformationObject.L

    # Compute number of dimensions:
    NumDimensions = length(X)

    # Generate samples of uncorrelated normal random variables U:
    USamples = [Distributions.rand(RNG, Distributions.Normal(), NumSamples, SamplingTechnique) for _ in 1:NumDimensions]
    USamples = vcat(USamples'...)

    # Generate samples of correlated normal random variables Z:
    ZSamples = L * USamples

    # Generate samples of correlated non-normal random variables X:
    XSamples = [Distributions.quantile.(X[i], Distributions.cdf.(Distributions.Normal(), ZSamples[i, :])) for i in 1:NumDimensions]
    XSamples = vcat(XSamples'...)

    return XSamples, ZSamples, USamples
end

Distributions.rand(TransformationObject::NatafTransformation, NumSamples::Integer, SamplingTechnique::Symbol) = 
Distributions.rand(Distributions.default_rng(), TransformationObject, NumSamples, SamplingTechnique)
