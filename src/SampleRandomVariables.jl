# Sample uncorrelated random variables:
function samplerv(Samplers::Union{<:Distribution,Vector{<:Distribution}}, NumSamples::Integer, SamplingTechnique::ITS)
    # Compute the number of distributions:
    NumDims = length(Samplers)

    # Preallocate:
    if NumDims == 1
        Samples = Vector{Float64}(undef, NumSamples)
    else
        Samples = Matrix{Float64}(undef, NumSamples, NumDims)
    end

    # Generate samples for each distribution:
    if NumDims == 1
        Samples = rand(Samplers, NumSamples)
    else
        for i = 1:NumDims
            Samples[:, i] = rand(Samplers[i], NumSamples)
        end
    end

    # Convert vector to scalar if only one sample if requested:
    if NumDims == 1 && NumSamples == 1
        Samples = Samples[1]
    end

    return Samples
end

function samplerv(Samplers::Union{<:Distribution,Vector{<:Distribution}}, NumSamples::Integer, SamplingTechnique::LHS)
    # Compute the number of distributions:
    NumDims = length(Samplers)

    # Preallocate:
    if NumDims == 1
        Samples = Vector{Float64}(undef, NumSamples)
    else
        Samples = Matrix{Float64}(undef, NumSamples, NumDims)
    end

    # Define the lower limits of each strata:
    LowerLimits = collect(range(0, (NumSamples - 1) / NumSamples, NumSamples))

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

    # Convert vector to scalar if only one sample if requested:
    if NumDims == 1 && NumSamples == 1
        Samples = Samples[1]
    end

    return Samples
end

# Sample correlated random variables:
function samplerv(Object::NatafTransformation, NumSamples::Integer)
    # Extract data:
    X = Object.X
    L = Object.L

    # Compute the number of marginal distributions:
    NumDims = length(X)

    # Generate samples of uncorrelated normal random variables U:
    USamples = Matrix{Float64}(undef, NumSamples, NumDims)
    for i in 1:NumDims
        USamples[:, i] = randn(NumSamples)
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