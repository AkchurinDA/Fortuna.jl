"""
    getdistortedcorrelation(X::Vector{<:Distribution}, ρˣ::Matrix{<:Real})

Returns a distorted correlation matrix ``\\underline{\\underline{\\rho}}^{Z}`` of correlated standard normal random variables ``\\underline{Z}``.
"""
function getdistortedcorrelation(X::Vector{<:Distribution}, ρˣ::Matrix{<:Real})
    # Compute the number of marginal distributions:
    NumDims = length(X)

    # Error-catching:
    if size(ρˣ) != (NumDims, NumDims)
        error("Size of the correlation matrix ρₓ is not compatible with the number of the marginal distributions.")
    end

    if !isposdef(ρˣ)
        error("Correlation matrix ρₓ must be a positive-definite matrix.")
    end

    MaxCorrelationValue = maximum(abs.(ρˣ - I))
    if MaxCorrelationValue > 1
        error("Off-diagonal entries of the correlation matrix ρˣ must be between -1 and +1.")
    end

    # Compute the locations and weights of the integration points in 1D:
    if MaxCorrelationValue < 0.9
        NumIP = 32
    else
        NumIP = 1024
    end

    IPLocations1D, IPWeights1D = gausslegendre(NumIP)

    # Transform the locations and weights of the integration points from 1D into 2D:
    ξ = Vector{Float64}(undef, NumIP^2)
    η = Vector{Float64}(undef, NumIP^2)
    Wᵢⱼ = Vector{Float64}(undef, NumIP^2)
    for i in 1:NumIP
        for j in 1:NumIP
            ξ[(i-1)*NumIP+j] = IPLocations1D[i]
            η[(i-1)*NumIP+j] = IPLocations1D[j]
            Wᵢⱼ[(i-1)*NumIP+j] = IPWeights1D[i] * IPWeights1D[j]
        end
    end

    # Set the bounds of integration:
    ZMin = -6
    ZMax = +6

    # Perform change of interval:
    zᵢ = ((ZMax - ZMin) / 2) * ξ .+ (ZMax + ZMin) / 2
    zⱼ = ((ZMax - ZMin) / 2) * η .+ (ZMax + ZMin) / 2

    # Compute the entries of the distorted correlation matrix:
    ρᶻ = Matrix{Float64}(I, NumDims, NumDims)
    for i in 1:NumDims
        for j in i+1:NumDims
            # Check if the marginal distributions are uncorrelated:
            if ρˣ[i, j] == 0
                continue
            end

            # Define a function from which we will compute the entries of the distorted correlation matrix:
            hᵢ = (quantile.(X[i], cdf.(Normal(0, 1), zᵢ)) .- mean(X[i])) / std(X[i])
            hⱼ = (quantile.(X[j], cdf.(Normal(0, 1), zⱼ)) .- mean(X[j])) / std(X[j])
            F(ρ, p) = ((ZMax - ZMin) / 2)^2 * dot(Wᵢⱼ .* (hᵢ .* hⱼ), ((1 / (2 * π * sqrt(1 - ρ^2))) * exp.((2 * ρ * (zᵢ .* zⱼ) - zᵢ .^ 2 - zⱼ .^ 2) / (2 * (1 - ρ^2))))) - ρˣ[i, j]

            # Compute the entries of the correlation matrix of the distorted correlation matrix:
            Problem = NonlinearProblem(F, ρˣ[i, j])
            Solution = solve(Problem, NewtonRaphson(), abstol=10^(-9), reltol=10^(-9))
            ρᶻ[i, j] = Solution.u
            ρᶻ[j, i] = ρᶻ[i, j]
        end
    end

    # Compute the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix and its inverse:
    C = cholesky(ρᶻ, check=true)
    L = C.L
    L⁻¹ = inv(L)

    return ρᶻ, L, L⁻¹
end

function transformsamples(Object::NatafTransformation, Samples::Union{Vector{<:Real},Matrix{<:Real}}, TransformationDirection::String)
    # Convert strings to lowercase:
    TransformationDirection = lowercase(TransformationDirection)

    if TransformationDirection != "x2u" && TransformationDirection != "u2x"
        error("Invalid transformation direction.")
    elseif TransformationDirection == "x2u"
        # Extract data:
        X = Object.X
        L⁻¹ = Object.L⁻¹

        # Check if the samples are passed as a vector or matrix:
        XSamples = Samples
        if typeof(XSamples) <: Vector
            XSamples = transpose(XSamples)
        end

        # Compute the number of samples and number of marginal distributions:
        NumSamples = size(XSamples, 1)
        NumDims = size(XSamples, 2)

        # Preallocate:
        ZSamples = Matrix{Float64}(undef, NumSamples, NumDims)

        # Convert samples of the marginal distributions X into the space of correlated standard normal random variables Z:
        for i in 1:NumDims
            ZSamples[:, i] = quantile.(Normal(0, 1), cdf.(X[i], XSamples[:, i]))
        end

        # Convert samples from the space of correlated standard normal random variables Z into the space of uncorrelated standard normal random variables U:
        USamples = transpose(L⁻¹ * transpose(ZSamples))

        # Check if the samples are passed as a vector or matrix:
        if typeof(transpose(XSamples)) <: Vector
            USamples = vec(USamples)
        end

        return USamples
    elseif TransformationDirection == "u2x"
        # Extract data:
        X = Object.X
        L = Object.L

        # Check if the samples are passed as a vector or matrix:
        USamples = Samples
        if typeof(USamples) <: Vector{<:Real}
            USamples = transpose(USamples)
        end

        # Compute the number of samples and number of marginal distributions:
        NumSamples = size(USamples, 1)
        NumDims = size(USamples, 2)

        # Convert samples to the space of correlated standard normal random variables Z:
        ZSamples = transpose(L * transpose(USamples))

        # Preallocate:
        XSamples = Matrix{Float64}(undef, NumSamples, NumDims)

        # Convert samples of the correlated standard normal random variables Z into samples of the marginals:
        for i in 1:NumDims
            XSamples[:, i] = quantile.(X[i], cdf.(Normal(0, 1), ZSamples[:, i]))
        end

        # Check if the samples are passed as a vector or matrix:
        if typeof(transpose(USamples)) <: Vector{<:Real}
            XSamples = vec(XSamples)
        end

        return XSamples
    end
end

function getjacobian(Object::NatafTransformation, Samples::Union{Vector{<:Real},Matrix{<:Real}}, TransformationDirection::String)
    # Convert strings to lowercase:
    TransformationDirection = lowercase(TransformationDirection)

    if TransformationDirection != "x2u" && TransformationDirection != "u2x"
        error("Invalid transformation direction.")
    elseif TransformationDirection == "x2u"
        # Extract data:
        X = Object.X
        L = Object.L

        # Check if the samples are passed as a vector or matrix:
        XSamples = Samples
        if typeof(XSamples) <: Vector
            # Convert column vector into row vector:
            XSamples = transpose(XSamples)
        end

        # Compute the number of samples and number of marginal distributions:
        NumSamples = size(XSamples, 1)
        NumDims = size(XSamples, 2)

        # Preallocate:
        ZSamples = Matrix{Float64}(undef, NumSamples, NumDims)

        # Convert samples to the space of correlated standard normal random variables Z:
        for i in 1:NumDims
            ZSamples[:, i] = quantile.(Normal(0, 1), cdf.(X[i], XSamples[:, i]))
        end

        # Preallocate:
        J = Array{Float64}(undef, NumDims, NumDims, NumSamples)

        # Compute the Jacobian:
        for i in 1:NumSamples
            D = zeros(NumDims, NumDims)
            for j in 1:NumDims
                D[j, j] = pdf(Normal(0, 1), ZSamples[i, j]) / pdf(X[j], XSamples[i, j])
            end
            J[:, :, i] = D * L
        end

        # Flatten out the 3D array into 2D array if the samples are passed as a vector:
        if NumSamples == 1
            J = J[:, :, 1]
        end
    elseif TransformationDirection == "u2x"
        # Extract data:
        X = Object.X
        L = Object.L
        L⁻¹ = Object.L⁻¹

        # Check if the samples are passed as a vector or matrix:
        USamples = Samples
        if typeof(USamples) <: Vector
            # Convert column vector into row vector:
            USamples = transpose(USamples)
        end

        # Compute the number of samples and number of marginal distributions:
        NumSamples = size(USamples, 1)
        NumDims = size(USamples, 2)

        # Convert samples to the space of correlated standard normal random variables Z:
        ZSamples = transpose(L * transpose(USamples))

        # Preallocate:
        XSamples = Matrix{Float64}(undef, NumSamples, NumDims)

        # Convert samples of the correlated standard normal random variables Z into samples of the marginals:
        for i in 1:NumDims
            XSamples[:, i] = quantile(X[i], cdf(Normal(0, 1), ZSamples[:, i]))
        end

        # Preallocate:
        J = Array{Float64}(undef, NumDims, NumDims, NumSamples)

        # Compute the Jacobian:
        for i in 1:NumSamples
            D = zeros(NumDims, NumDims)
            for j in 1:NumDims
                D[j, j] = pdf(X[j], XSamples[i, j]) / pdf(Normal(0, 1), ZSamples[i, j])
            end
            J[:, :, i] = L⁻¹ * D
        end

        # Flatten out the 3D array into 2D array if the samples are passed as a vector:
        if NumSamples == 1
            J = J[:, :, 1]
        end
    end

    return J
end

"""
    jointpdf(Object::NatafTransformation, XSamples::Union{Vector{<:Real},Matrix{<:Real}})

Returns values of joint probability density functions ``f_{\\underline{X}}(\\undeline{x})`` of non-normal correlated random variables ``\\underline{X}`` evaluated at given points.
"""
function jointpdf(Object::NatafTransformation, XSamples::Union{Vector{<:Real},Matrix{<:Real}})
    # Extract data:
    X = Object.X
    ρᶻ = Object.ρᶻ

    # Check if the samples are passed as a vector or matrix:
    if typeof(XSamples) <: Vector
        # Convert column vector into row vector:
        XSamples = transpose(XSamples)
    end

    # Compute the number of samples and number of marginal distributions:
    NumSamples = size(XSamples, 1)
    NumDims = size(XSamples, 2)

    # Preallocate:
    ZSamples = Matrix{Float64}(undef, NumSamples, NumDims)
    f = Matrix{Float64}(undef, NumSamples, NumDims)
    ϕ = Matrix{Float64}(undef, NumSamples, NumDims)

    # Convert samples to the space of correlated standard normal random variables Z:
    for i in 1:NumDims
        ZSamples[:, i] = quantile(Normal(0, 1), cdf(X[i], XSamples[:, i]))
        f[:, i] = pdf.(X[i], XSamples[:, i])
        ϕ[:, i] = pdf.(Normal(0, 1), ZSamples[:, i])
    end

    # Compute the joint PDF of samples in the space of correlated standard normal random variables Z: 
    JointPDFZ = pdf(MvNormal(ρᶻ), transpose(ZSamples))

    # Compute the joint PDF:
    JointPDFX = JointPDFZ .* (prod(f, dims=2) ./ prod(ϕ, dims=2))

    # Clean up:
    JointPDFX = replace(JointPDFX, NaN => 0, Inf => 0)
    JointPDFX = vec(JointPDFX)

    return JointPDFX
end