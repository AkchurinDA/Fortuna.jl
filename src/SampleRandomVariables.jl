# Sample random vectors with uncorrelated marginal random variables:
"""
    samplerv(Samplers::Union{<:Distribution,Vector{<:Distribution}}, NumSamples::Integer, SamplingTechnique::AbstractSamplingTechnique)

The function returns samples of random variables and random vectors with uncorrelated marginal random variables using various sampling techniques.
- If `SamplingTechnique = ITS()`, the function generates samples using Inverse Transform Sampling technique.
- If `SamplingTechnique = LHS()`, the function generates samples using Latin Hypercube Sampling technique.
"""
function samplerv(Samplers::Union{<:Distribution, Vector{<:Distribution}}, NumSamples::Integer, SamplingTechnique::AbstractSamplingTechnique)
    # Compute the number of distributions:
    NumDims = length(Samplers)

    # Preallocate:
    if NumDims == 1
        Samples = Vector{Float64}(undef, NumSamples)
    else
        Samples = Matrix{Float64}(undef, NumSamples, NumDims)
    end

    if isa(SamplingTechnique, ITS)
        # Generate samples for each distribution:
        if NumDims == 1
            Samples = rand(Samplers, NumSamples)
        else
            for i = 1:NumDims
                Samples[:, i] = rand(Samplers[i], NumSamples)
            end
        end
    elseif isa(SamplingTechnique, LHS)
        # Define the lower limits of each strata:
        LowerLimits = collect(0:(1/NumSamples):((NumSamples-1)/NumSamples))

        if NumDims == 1
            # Generate samples from a uniform distributions:
            UniformSamples = LowerLimits + rand(Uniform(0, 1 / NumSamples), NumSamples)

            # Shuffle samples from a uniform distributions:
            UniformSamples = shuffle(UniformSamples)

            # Generate samples:
            Samples = quantile.(Samplers, UniformSamples)
        else
            for i = 1:NumDims
                # Generate samples from a uniform distributions:
                UniformSamples = LowerLimits + rand(Uniform(0, 1 / NumSamples), NumSamples)

                # Shuffle samples from a uniform distributions:
                UniformSamples = shuffle(UniformSamples)

                # Generate samples:
                Samples[:, i] = quantile.(Samplers[i], UniformSamples)
            end
        end
    else
        error("Provided sampling technique is not supported.")
    end

    # Convert vector to scalar if only one sample if requested:
    if NumDims == 1 && NumSamples == 1
        Samples = Samples[1]
    end

    return Samples
end

# Sample random vectors with correlated marginal random variables:
"""
    samplerv(Object::NatafTransformation, NumSamples::Integer, SamplingTechnique::AbstractSamplingTechnique)

This function generates samples of random variables in ``X``-, ``Z``-, and ``U``-spaces using a `NatafTransformation` object.
- ``X``-space - space of correlated non-normal random variables
- ``Z``-space - space of correlated standard normal random variables
- ``U``-space - space of uncorrelated standard normal random variables
"""
function samplerv(Object::NatafTransformation, NumSamples::Integer, SamplingTechnique::AbstractSamplingTechnique)
    # Extract data:
    X = Object.X
    L = Object.L

    # Compute the number of marginal distributions:
    NumDims = length(X)

    # Define a standard normal random variable:
    U = generaterv("Normal", "M", [0, 1])

    # Generate samples of uncorrelated normal random variables U:
    USamples = Matrix{Float64}(undef, NumSamples, NumDims)
    for i in 1:NumDims
        USamples[:, i] = samplerv(U, NumSamples, SamplingTechnique)
    end

    # Generate samples of correlated normal random variables Z:
    ZSamples = USamples * transpose(L)

    # Generate samples of correlated non-normal random variables X:
    XSamples = Matrix{Float64}(undef, NumSamples, NumDims)
    for i in 1:NumDims
        XSamples[:, i] = quantile.(X[i], cdf.(Normal(0, 1), ZSamples[:, i]))
    end

    return XSamples, ZSamples, USamples
end