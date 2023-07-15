# Sample uncorrelated random variables:
function samplerv(Samplers::Union{Distribution,Vector{<:Distribution}}, NumSamples::Integer, SamplingTechnique::String)
    # Convert strings to lowercase:
    SamplingTechnique = lowercase(SamplingTechnique)

    # Compute the number of distributions:
    NumDims = length(Samplers)

    # Preallocate:
    if NumDims == 1
        Samples = Vector{Float64}(undef, NumSamples)
    else
        Samples = Matrix{Float64}(undef, NumSamples, NumDims)
    end

    if SamplingTechnique == "inversetransformsampling" || SamplingTechnique == "its"
        if NumDims == 1
            # Generate samples from a uniform distributions:
            UniformSamples = rand(NumSamples)

            # Generate samples for each distribution:
            Samples = quantile(Samplers, UniformSamples)
        else
            for i = 1:NumDims
                # Generate samples from a uniform distributions:
                UniformSamples = rand(NumSamples)

                # Generate samples for each distribution:
                Samples[:, i] = quantile(Samplers[i], UniformSamples)
            end
        end
    elseif SamplingTechnique == "latinhypercubesampling" || SamplingTechnique == "lhs"
        # Define the lower limits of each strata:
        LowerLimits = collect(range(0, (NumSamples - 1) / NumSamples, NumSamples))

        if NumDims == 1
            # Generate samples from a uniform distributions:
            UniformSamples = LowerLimits + rand(Uniform(0, 1 / NumSamples), NumSamples)

            # Shuffle samples from a uniform distributions:
            UniformSamples = shuffle(UniformSamples)

            # Generate samples:
            Samples = quantile(Samplers, UniformSamples)
        else
            for i = 1:NumDims
                # Generate samples from a uniform distributions:
                UniformSamples = LowerLimits + rand(Uniform(0, 1 / NumSamples), NumSamples)

                # Shuffle samples from a uniform distributions:
                UniformSamples = shuffle(UniformSamples)

                # Generate samples:
                Samples[:, i] = quantile(Samplers[i], UniformSamples)
            end
        end
    else
        error("Sampling technique is not defined.")
    end

    # Convert vector to scalar if only one sample if requested:
    if NumDims == 1 && NumSamples == 1
        Samples = Samples[1]
    end

    return Samples
end

# Sample correlated random variables:
function samplerv(Object::NatafTransformation, NumSamples::Integer, SamplingTechnique::String)
    # Extract data:
    X = Object.X
    L = Object.L

    # Compute the number of marginal distributions:
    NumDims = length(X)

    # Generate samples of uncorrelated normal random variables U:
    NormalDistribution = generaterv("Normal", "Moments", [0, 1])
    USamples = Matrix{Float64}(undef, NumSamples, NumDims)
    for i in 1:NumDims
        USamples[:, i] = samplerv(NormalDistribution, NumSamples, SamplingTechnique)
    end

    # Generate samples of correlated normal random variables Z:
    ZSamples = USamples * transpose(L)

    # Generate samples of correlated non-normal random variables X:
    XSamples = Matrix{Float64}(undef, NumSamples, NumDims)
    for i in 1:NumDims
        XSamples[:, i] = quantile(X[i], cdf(Normal(0, 1), ZSamples[:, i]))
    end

    return XSamples, ZSamples, USamples
end