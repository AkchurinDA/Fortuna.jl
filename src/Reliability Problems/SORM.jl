"""
    SORM <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using Second-Order Reliability Method (SORM).

$(TYPEDFIELDS)
"""
Base.@kwdef struct SORM <: AbstractReliabililyAnalysisMethod
    Submethod::SORMSubmethod = CF()
end

"""
    CF <: SORMSubmethod

Type used to perform reliability analysis using Curve-Fitting (CF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct CF <: SORMSubmethod # Curve-Fitting method
    "Step size used to compute the Hessian at the design point in ``U``-space"
    ϵ::Real = 1 / 1000
end

"""
    PF <: SORMSubmethod

Type used to perform reliability analysis using Point-Fitting (PF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct PF <: SORMSubmethod # Point-Fitting method

end

"""
    CFCache

Type used to perform reliability analysis using Point-Fitting (PF) method.

$(TYPEDFIELDS)
"""
struct CFCache # Curve-Fitting method
    "Results of reliability analysis performed using First-Order Reliability Method (FORM)"
    FORMSolution    ::iHLRFCache
    "Generalized reliability indices ``\\beta``"
    β₂              ::Vector{Float64}
    "Probabilities of failure ``P_{f}``"
    PoF₂            ::Vector{Float64}
    "Principal curvatures ``\\kappa``"
    κ               ::Vector{Float64}
end

"""
    PFCache

Type used to perform reliability analysis using Point-Fitting (PF) method.

$(TYPEDFIELDS)
"""
struct PFCache # Point-Fitting method
    "Results of reliability analysis performed using First-Order Reliability Method (FORM)"
    FORMSolution    ::iHLRFCache
    "Generalized reliability index ``\\beta``"
    β₂              ::Vector{Float64}
    "Probabilities of failure ``P_{f}``"
    PoF₂            ::Vector{Float64}
    "Fitting points on the negative side of the hyper-cylinder"
    FittingPoints⁻  ::Matrix{Float64}
    "Fitting points on the positive side of the hyper-cylinder"
    FittingPoints⁺  ::Matrix{Float64}
    "Principal curvatures on the negative and positive sides"
    κ₁              ::Matrix{Float64}
    "Principal curvatures of each hyper-semiparabola"
    κ₂              ::Matrix{Float64}
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::SORM)

Function used to solve reliability problems using Second-Order Reliability Method (SORM).
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::SORM)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    # Determine the design point using FORM:
    FORMSolution    = solve(Problem, FORM())
    u               = FORMSolution.u[:, end]
    ∇G              = FORMSolution.∇G[:, end]
    α               = FORMSolution.α[:, end]
    β₁              = FORMSolution.β

    # Extract the problem data:
    X   = Problem.X
    ρˣ  = Problem.ρˣ
    g   = Problem.g

    if !isa(Submethod, CF) && !isa(Submethod, PF)
        error("Invalid SORM submethod.")
    elseif isa(Submethod, CF)
        # Extract the analysis details:
        ϵ = Submethod.ϵ

        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Perform Nataf transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Compute the Hessian at the design point in U-space:
        H = gethessian(g, NatafObject, NumDimensions, u, ϵ)

        # Compute the orthonomal matrix:
        R = getorthonormal(α, NumDimensions)

        # Evaluate the principal curvatures:
        A = R * H * LinearAlgebra.transpose(R) / LinearAlgebra.norm(∇G)
        κ = LinearAlgebra.eigen(A[1:end-1, 1:end-1]).values

        # Compute the probabilities of failure:
        PoF₂ = Vector{Float64}(undef, 2)

        begin # Hohenbichler-Rackwitz (1988)
            ψ = Distributions.pdf(Distributions.Normal(), β₁) / Distributions.cdf(Distributions.Normal(), -β₁)

            if all(κᵢ -> ψ * κᵢ > -1, κ)
                PoF₂[1] = Distributions.cdf(Distributions.Normal(), -β₁) * prod(κᵢ -> 1 / sqrt(1 + ψ * κᵢ), κ)
            else
                PoF₂[1] = nothing
                error("Condition of Hohenbichler-Rackwitz's approximation of the probability of failure was not satisfied.")
            end
        end

        begin # Breitung (1984)
            if all(κᵢ -> β₁ * κᵢ > -1, κ)
                PoF₂[2] = Distributions.cdf(Distributions.Normal(), -β₁) * prod(κᵢ -> 1 / sqrt(1 + β₁ * κᵢ), κ)
            else
                PoF₂[2] = nothing
                error("Condition of Breitung's approximation of the probability of failure was not satisfied.")
            end
        end

        # Compute the generalized reliability index:
        β₂ = -Distributions.quantile.(Distributions.Normal(), PoF₂)

        # Return results:
        return CFCache(FORMSolution, β₂, PoF₂, κ)
    elseif isa(Submethod, PF)
        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Perform Nataf transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Compute the orthonomal matrix:
        R = getorthonormal(α, NumDimensions)
        u′ = R * u
        if u′[end] < 0
            R = -R
        end

        # Compute radius of a hypercylinder:
        if β₁ < 1
            H = 1
        elseif β₁ ≥ 1 && β₁ ≤ 3
            H = β₁
        elseif β₁ > 3
            H = 3
        end

        # Compute fitting points:
        FittingPoints⁺  = Matrix{Float64}(undef, NumDimensions - 1, 2)
        FittingPoints⁻  = Matrix{Float64}(undef, NumDimensions - 1, 2)
        κ₁              = Matrix{Float64}(undef, NumDimensions - 1, 2)
        for i in 1:(NumDimensions - 1)
            function F(u, p)
                UPrime          = zeros(eltype(u), NumDimensions)
                UPrime[i]       = p
                UPrime[end]     = u
            
                return G′(g, NatafObject, R, UPrime)
            end

            # Negative side:
            Problem⁻                = NonlinearSolve.NonlinearProblem(F, β₁, -H)
            Solution⁻               = NonlinearSolve.solve(Problem⁻, NonlinearSolve.NewtonRaphson(), abstol=10^(-9), reltol=10^(-9))
            FittingPoints⁻[i, 1]    = -H
            FittingPoints⁻[i, 2]    = Solution⁻.u

            # Positive side:
            Problem⁺                = NonlinearSolve.NonlinearProblem(F, β₁, +H)
            Solution⁺               = NonlinearSolve.solve(Problem⁺, NonlinearSolve.NewtonRaphson(), abstol=10^(-9), reltol=10^(-9))
            FittingPoints⁺[i, 1]    = +H
            FittingPoints⁺[i, 2]    = Solution⁺.u

            # Curvatures:
            κ₁[i, 1] = 2 * (FittingPoints⁻[i, 2] - β₁) / (FittingPoints⁻[i, 1] ^ 2) # Negative side
            κ₁[i, 2] = 2 * (FittingPoints⁺[i, 2] - β₁) / (FittingPoints⁺[i, 1] ^ 2) # Positive side
        end

        # Compute number of hyperquadrants used to fit semiparabolas:
        NumHyperquadrants = 2 ^ (NumDimensions - 1)

        # Get all possible permutations of curvatures:
        Indices = Base.Iterators.repeated(1:2, NumDimensions - 1)
        Indices = Base.Iterators.product(Indices...)
        Indices = collect(Indices)
        Indices = vec(Indices)

        κ₂ = Matrix{Float64}(undef, NumHyperquadrants, NumDimensions - 1)
        for i in 1:NumHyperquadrants
            for j in 1:(NumDimensions - 1)
                κ₂[i, j] = κ₁[j, Indices[i][j]]
            end
        end

        # Compute the probabilities of failure for each semiparabola:
        PoF₂ = Matrix{Float64}(undef, NumHyperquadrants, 2)
        for i in 1:NumHyperquadrants
            κ = κ₂[i, :]

            begin # Hohenbichler-Rackwitz (1988)
                ψ = Distributions.pdf(Distributions.Normal(), β₁) / Distributions.cdf(Distributions.Normal(), -β₁)

                if all(κᵢ -> ψ * κᵢ > -1, κ)
                    PoF₂[i, 1] = Distributions.cdf(Distributions.Normal(), -β₁) * prod(κᵢ -> 1 / sqrt(1 + ψ * κᵢ), κ)
                else
                    error("Condition of Hohenbichler-Rackwitz's approximation of the probability of failure was not satisfied.")
                end
            end

            begin # Breitung (1984)
                if all(κᵢ -> β₁ * κᵢ > -1, κ)
                    PoF₂[i, 2] = Distributions.cdf(Distributions.Normal(), -β₁) * prod(κᵢ -> 1 / sqrt(1 + β₁ * κᵢ), κ)
                else
                    error("Condition of Breitung's approximation of the probability of failure was not satisfied.")
                end
            end
        end

        PoF₂ = (1/ NumHyperquadrants) * PoF₂
        PoF₂ = sum(PoF₂, dims = 1)
        PoF₂ = vec(PoF₂)

        # Compute the generalized reliability index:
        β₂ = -Distributions.quantile.(Distributions.Normal(), PoF₂)

        # Return results:
        return PFCache(FORMSolution, β₂, PoF₂, FittingPoints⁻, FittingPoints⁺, κ₁, κ₂)
    end
end

function gethessian(g::Function, NatafObject::NatafTransformation, NumDimensions::Integer, u::Vector{Float64}, ϵ::Real)
    # Preallocate:
    H = Matrix{Float64}(undef, NumDimensions, NumDimensions)

    for i in 1:NumDimensions
        for j in 1:NumDimensions
            # Define the pertubation directions:
            eᵢ = zeros(NumDimensions,)
            eⱼ = zeros(NumDimensions,)
            eᵢ[i] = 1
            eⱼ[j] = 1

            # Perturb the design point in U-space:
            u₁ = u + ϵ * eᵢ + ϵ * eⱼ
            u₂ = u + ϵ * eᵢ - ϵ * eⱼ
            u₃ = u - ϵ * eᵢ + ϵ * eⱼ
            u₄ = u - ϵ * eᵢ - ϵ * eⱼ

            # Transform the perturbed design points from X- to U-space:
            x₁ = transformsamples(NatafObject, u₁, "U2X")
            x₂ = transformsamples(NatafObject, u₂, "U2X")
            x₃ = transformsamples(NatafObject, u₃, "U2X")
            x₄ = transformsamples(NatafObject, u₄, "U2X")

            # Evaluate the limit state function at the perturbed points:
            G₁ = g(x₁)
            G₂ = g(x₂)
            G₃ = g(x₃)
            G₄ = g(x₄)

            # Evaluate the entries of the Hessian using finite difference method:
            H[i, j] = (G₁ - G₂ - G₃ + G₄) / (4 * ϵ^2)
        end
    end

    return H
end

function getorthonormal(α::Vector{Float64}, NumDimensions::Integer)
    # Initilize the matrix:
    A = Matrix(1.0 * I, NumDimensions, NumDimensions)
    A = reverse(A, dims=2)
    A[:, 1] = LinearAlgebra.transpose(α)

    # Perform QR factorization:
    Q, _ = LinearAlgebra.qr(A)
    Q = Matrix(Q)

    # Clean up the result:
    R = LinearAlgebra.transpose(reverse(Q, dims=2))
    R = Matrix(R)

    return R
end

function G′(g::Function, NatafObject::NatafTransformation, R::Matrix{Float64}, UPrimeSamples::AbstractVector)
    # Transform samples from U'- to X-space:
    USamples = LinearAlgebra.transpose(R) * UPrimeSamples
    XSamples = transformsamples(NatafObject, USamples, "U2X")

    # Evaluate the limit state function at the transform samples:
    GPrimeSamples = g(XSamples)

    # Return the result:
    return GPrimeSamples
end