"""
    NatafTransformation <: AbstractTransformation

Type used to perform Nataf Transformation.

$(TYPEDFIELDS)
"""
mutable struct NatafTransformation <: AbstractIsoprobabilisticTransformation
    "Random vector ``\\vec{X}``"
    X   ::AbstractVector{<:Distributions.ContinuousUnivariateDistribution}
    "Correlation matrix ``\\rho^{X}``"
    ρˣ  ::AbstractMatrix{<:Real}
    "Distorted correlation matrix ``\\rho^{Z}``"
    ρᶻ  ::AbstractMatrix{<:Real}
    "Lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix ``L``"
    L   ::AbstractMatrix{<:Real}
    "Inverse of the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix ``L^{-1}``"
    L⁻¹ ::AbstractMatrix{<:Real}

    function NatafTransformation(X::AbstractVector{<:Distributions.ContinuousUnivariateDistribution}, ρˣ::AbstractMatrix{<:Real})
        # Compute the distorted correlation matrix:
        ρᶻ, L, L⁻¹ = getdistortedcorrelation(X, ρˣ)

        # Return the Nataf Transformation object with the computed distorted correlation matrix:
        return new(X, ρˣ, ρᶻ, L, L⁻¹)
    end
end
Base.broadcastable(x::NatafTransformation) = Ref(x)

"""
    getdistortedcorrelation(X::AbstractVector{<:Distributions.UnivariateDistribution}, ρˣ::AbstractMatrix{<:Real})

Function used to compute the distorted correlation matrix ``\\rho^{Z}``.
"""
function getdistortedcorrelation(X::AbstractVector{<:Distributions.ContinuousUnivariateDistribution}, ρˣ::AbstractMatrix{<:Real}) 
    # Compute number of dimensions:
    NumDimensions = length(X)

    # Error-catching:
    if size(ρˣ) != (NumDimensions, NumDimensions)
        error("Size of the correlation matrix ρₓ is not compatible with the number of marginal distributions.")
    end

    if !LinearAlgebra.isposdef(ρˣ)
        error("Correlation matrix ρₓ must be a positive-definite matrix.")
    end

    MaxCorrelationValue = maximum(abs.(ρˣ - I))
    if MaxCorrelationValue > 1
        error("Off-diagonal entries of the correlation matrix ρˣ must be between -1 and +1.")
    end

    # Compute the locations and weights of the integration points in 1D:
    NumIP                      = MaxCorrelationValue ≤ 0.9 ? 64 : 1024
    IPLocations1D, IPWeights1D = FastGaussQuadrature.gausslegendre(NumIP)

    # Transform the locations and weights of the integration points from 1D into 2D:
    ξ = Vector{Float64}(undef, NumIP ^ 2)
    η = Vector{Float64}(undef, NumIP ^ 2)
    W = Vector{Float64}(undef, NumIP ^ 2)
    for i in 1:NumIP
        for j in 1:NumIP
            ξ[(i - 1) * NumIP + j] = IPLocations1D[i]
            η[(i - 1) * NumIP + j] = IPLocations1D[j]
            W[(i - 1) * NumIP + j] = IPWeights1D[i] * IPWeights1D[j]
        end
    end

    # Set the bounds of integration:
    ZMin = -6
    ZMax = +6

    # Perform change of interval:
    zᵢ = ((ZMax - ZMin) / 2) * ξ .+ (ZMax + ZMin) / 2
    zⱼ = ((ZMax - ZMin) / 2) * η .+ (ZMax + ZMin) / 2

    # Determine the common parameter type:
    Parameters          = params.(X)
    Parameters          = [Parameters[i][j] for i in eachindex(X) for j in eachindex(Parameters[i])]
    ParameterTypes      = typeof.(Parameters)
    CommonParameterType = promote_type(ParameterTypes...)

    # Compute the entries of the distorted correlation matrix:
    ρᶻ = Matrix{CommonParameterType}(I, NumDimensions, NumDimensions)
    for i in 1:NumDimensions
        for j in (i + 1):NumDimensions
            # Check if the marginal distributions are uncorrelated:
            if ρˣ[i, j] == 0
                continue
            end

            # Define a function from which we will compute the entries of the distorted correlation matrix:
            hᵢ      = (Distributions.quantile.(X[i], Distributions.cdf.(Distributions.Normal(), zᵢ)) .- Distributions.mean(X[i])) / Distributions.std(X[i])
            hⱼ      = (Distributions.quantile.(X[j], Distributions.cdf.(Distributions.Normal(), zⱼ)) .- Distributions.mean(X[j])) / Distributions.std(X[j])
            F(ρ, p) = ((ZMax - ZMin) / 2) ^ 2 * LinearAlgebra.dot(W .* (hᵢ .* hⱼ), ((1 / (2 * π * sqrt(1 - ρ ^ 2))) * exp.((2 * ρ * (zᵢ .* zⱼ) - zᵢ .^ 2 - zⱼ .^ 2) / (2 * (1 - ρ ^ 2))))) - ρˣ[i, j]

            # Compute the entries of the correlation matrix of the distorted correlation matrix:
            try
                Problem  = NonlinearSolve.NonlinearProblem(F, ρˣ[i, j])
                Solution = NonlinearSolve.solve(Problem, nothing)
                ρᶻ[i, j] = Solution.u
            catch
                Problem  = NonlinearSolve.IntervalNonlinearProblem(F, (-(1 - 1E-3), +(1 - 1E-3)))
                Solution = NonlinearSolve.solve(Problem, nothing)
                ρᶻ[i, j] = Solution.u
            end
            ρᶻ[j, i] = ρᶻ[i, j]
        end
    end

    # Compute the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix and its inverse:
    L   = LinearAlgebra.cholesky(ρᶻ).L
    L⁻¹ = LinearAlgebra.inv(L)

    # Return the result:
    return ρᶻ, L, L⁻¹
end

"""
    transformsamples(TransformationObject::NatafTransformation, Samples::AbstractVector{<:Real}, TransformationDirection::Symbol)

Function used to transform samples from ``X``- to ``U``-space and vice versa. \\
If `TransformationDirection is:
- `:X2U`, then the function transforms samples ``\\vec{x}`` from ``X``- to ``U``-space.
- `:U2X`, then the function transforms samples ``\\vec{u}`` from ``U``- to ``X``-space.
"""
function transformsamples(TransformationObject::NatafTransformation, Samples::AbstractVector{<:Real}, TransformationDirection::Symbol)
    # Compute number of dimensions:
    NumDimensions = length(Samples)

    if TransformationDirection != :X2U && TransformationDirection != :U2X
        error("Invalid transformation direction! Available options are: :X2U and :U2X.")
    elseif TransformationDirection == :X2U
        # Extract data:
        X   = TransformationObject.X
        L⁻¹ = TransformationObject.L⁻¹

        # Convert samples of the marginal distributions X into the space of correlated standard normal random variables Z:
        ZSamples = [Distributions.quantile(Distributions.Normal(), Distributions.cdf(X[i], Samples[i])) for i in 1:NumDimensions]
        
        # Convert samples from the space of correlated standard normal random variables Z into the space of uncorrelated standard normal random variables U:
        USamples = L⁻¹ * ZSamples

        # Return the result:
        return USamples
    elseif TransformationDirection == :U2X
        # Extract data:
        X = TransformationObject.X
        L = TransformationObject.L

        # Convert samples to the space of correlated standard normal random variables Z:
        ZSamples = L * Samples

        # Convert samples of the correlated standard normal random variables Z into samples of the marginals:
        XSamples = [Distributions.quantile(X[i], Distributions.cdf(Distributions.Normal(), ZSamples[i])) for i in 1:NumDimensions]

        # Return the result:
        return XSamples
    end
end

function transformsamples(TransformationObject::NatafTransformation, Samples::AbstractMatrix{<:Real}, TransformationDirection::Symbol)
    # Compute number of samples and dimensions:
    NumDimensions = size(Samples, 1)
    NumSamples    = size(Samples, 2)

    # Error-catching:
    if NumDimensions != length(TransformationObject.X)
        error("Number of dimensions of the samples is incorrect!")
    end

    TransformedSamples = [transformsamples(TransformationObject, Samples[:, i], TransformationDirection) for i in 1:NumSamples]
    TransformedSamples = hcat(TransformedSamples...)

    # Return the result:
    return TransformedSamples
end

"""
    getjacobian(TransformationObject::NatafTransformation, Samples::AbstractVector{<:Real}, TransformationDirection::Symbol)

Function used to compute the Jacobians of the transformations of samples from ``X``- to ``U``-space and vice versa. \\
If `TransformationDirection` is:
- `:X2U`, then the function returns the Jacobians of the transformations of samples ``\\vec{x}`` from ``X``- to ``U``-space.
- `:U2X`, then the function returns the Jacobians of the transformations of samples ``\\vec{u}`` from ``U``- to ``X``-space.
"""
function getjacobian(TransformationObject::NatafTransformation, Samples::AbstractVector{<:Real}, TransformationDirection::Symbol)
    # # Convert strings to lowercase:
    # TransformationDirection = lowercase(TransformationDirection)

    # Compute number of dimensions:
    NumDimensions = length(Samples)

    if TransformationDirection != :X2U && TransformationDirection != :U2X
        error("Invalid transformation direction.")
    elseif TransformationDirection == :X2U
        # Extract data:
        X = TransformationObject.X
        L = TransformationObject.L

        # Convert samples to the space of correlated standard normal random variables Z:
        ZSamples = [Distributions.quantile(Distributions.Normal(), Distributions.cdf(X[i], Samples[i])) for i in 1:NumDimensions]

        # Compute the Jacobian:
        D = [Distributions.pdf(Distributions.Normal(), ZSamples[i]) / Distributions.pdf(X[i], Samples[i]) for i in 1:NumDimensions]
        D = LinearAlgebra.diagm(D)
        J = D * L
    elseif TransformationDirection == :U2X
        # Extract data:
        X   = TransformationObject.X
        L   = TransformationObject.L
        L⁻¹ = TransformationObject.L⁻¹

        # Convert samples to the space of correlated standard normal random variables Z:
        ZSamples = L * Samples

        # Convert samples of the correlated standard normal random variables Z into samples of the marginals:
        XSamples = [Distributions.quantile(X[i], Distributions.cdf(Distributions.Normal(), ZSamples[i])) for i in 1:NumDimensions]

        # Compute the Jacobian:
        D = [Distributions.pdf(X[i], XSamples[i]) / Distributions.pdf(Distributions.Normal(), ZSamples[i]) for i in 1:NumDimensions]
        D = LinearAlgebra.diagm(D)
        J = L⁻¹ * D
    end

    # Return the result:
    return J
end

function getjacobian(TransformationObject::NatafTransformation, Samples::AbstractMatrix{<:Real}, TransformationDirection::Symbol)
    # Compute number of dimensions and samples:
    NumDimensions = size(Samples, 1)
    NumSamples    = size(Samples, 2)

    # Error-catching:
    if NumDimensions != length(TransformationObject.X)
        error("Number of dimensions of the samples is incorrect.")
    end

    # Compute the Jacobians:
    J = [getjacobian(TransformationObject, Samples[:, i], TransformationDirection) for i in 1:NumSamples]

    # Return the result:
    return J
end

"""
    pdf(TransformationObject::NatafTransformation, x::AbstractVector{<:Real})

Function used to compute the joint PDF in ``X``-space.
"""
function Distributions.pdf(TransformationObject::NatafTransformation, x::AbstractVector{<:Real})
    # Extract data:
    X  = TransformationObject.X
    ρᶻ = TransformationObject.ρᶻ

    # Compute the number of samples and number of marginal distributions:
    NumDimensions = length(x)

    # Convert samples to the space of correlated standard normal random variables Z:
    z = [Distributions.quantile(Distributions.Normal(), Distributions.cdf(X[i], x[i])) for i in 1:NumDimensions]
    f = [Distributions.pdf(X[i], x[i])                                                 for i in 1:NumDimensions]
    ϕ = [Distributions.pdf(Distributions.Normal(), z[i])                               for i in 1:NumDimensions]

    # Compute the joint PDF of samples in the space of correlated standard normal random variables Z: 
    JointPDFZ = Distributions.pdf(Distributions.MvNormal(ρᶻ), z)

    # Compute the joint PDF:
    JointPDFX = JointPDFZ * (prod(f) / prod(ϕ))

    # Clean up:
    if !LinearAlgebra.isfinite(JointPDFX)
        JointPDFX = 0
    end

    # Return the result:
    return JointPDFX
end

function Distributions.pdf(TransformationObject::NatafTransformation, x::AbstractMatrix{<:Real})
    # Compute number of dimensions and samples:
    NumDimensions = size(x, 1)
    NumSamples    = size(x, 2)

    # Error-catching:
    if NumDimensions != length(TransformationObject.X)
        error("Number of dimensions of the samples is incorrect.")
    end

    # Compute the joint PDF:
    JointPDFX = [Distributions.pdf(TransformationObject, x[:, i]) for i in 1:NumSamples]

    # Return the result:
    return JointPDFX
end