# Second-Order Reliability Method:
"""
    analyze(Problem::ReliabilityProblem, AnalysisMethod::SORM)

The function solves the provided reliability problem using any submethod that falls under a broader category of Second-Order Reliability Methods (SORM).
"""
function analyze(Problem::ReliabilityProblem, AnalysisMethod::SORM)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    # Determine the design point using FORM:
    FORMSolution = analyze(Problem, FORM(iHLRF()))
    β₁ = FORMSolution.β
    PoF₁ = FORMSolution.PoF
    u = FORMSolution.u[:, end]
    ∇G = FORMSolution.∇G[:, end]
    α = FORMSolution.α[:, end]

    # Extract the problem data:
    g = Problem.g
    X = Problem.X
    ρˣ = Problem.ρˣ

    if !isa(Submethod, CF) && !isa(Submethod, GF) && !isa(Submethod, PF)
        error("Invalid SORM submethod.")
    elseif isa(Submethod, CF)
        # Extract the analysis details:
        ϵ = Submethod.ϵ

        # Compute the number of marginal distributions:
        NumDims = length(X)

        # Perform Nataf transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Compute the Hessian at the design point in U-space:
        H = gethessian(NatafObject, g, u, ϵ, NumDims)

        # Compute the orthonomal matrix:
        R = getorthonormal(α, NumDims)

        # Evaluate the principal curvatures:
        A = R * H * transpose(R) / norm(∇G)
        κ = eigen(A[1:end-1, 1:end-1]).values

        # Compute the probabilities of failure:
        PoF₂ = Vector{Float64}(undef, 2)

        begin # Hohenbichler-Rackwitz (1988)
            ψ = pdf(Normal(0, 1), β₁) / cdf(Normal(0, 1), -β₁)

            if all(κᵢ -> ψ * κᵢ > -1, κ)
                PoF₂[1] = cdf(Normal(0, 1), -β₁) * prod(κᵢ -> 1 / sqrt(1 + ψ * κᵢ), κ)
            else
                PoF₂[1] = nothing
                println("Condition of Hohenbichler-Rackwitz's approximation of the probability of failure was not satisfied.")
            end
        end

        begin # Breitung (1984)
            if all(κᵢ -> β₁ * κᵢ > -1, κ)
                PoF₂[2] = cdf(Normal(0, 1), -β₁) * prod(κᵢ -> 1 / sqrt(1 + β₁ * κᵢ), κ)
            else
                PoF₂[2] = nothing
                println("Condition of Breitung's approximation of the probability of failure was not satisfied.")
            end
        end

        # Compute the generalized reliability index:
        β₂ = -quantile.(Normal(0, 1), PoF₂)

        # Return results:
        return CFCache(β₁, PoF₁, β₂, PoF₂, H, R, A, κ)
    elseif isa(Submethod, GF)
        # Not yet implemented
    elseif isa(Submethod, PF)
        # Not yet implemented
    end
end

function gethessian(NatafObject, g, u, ϵ, NumDims)
    # Preallocate:
    H = Matrix{Float64}(undef, NumDims, NumDims)

    for i in 1:NumDims
        for j in 1:NumDims
            # Define the pertubation directions:
            eᵢ = zeros(NumDims,)
            eⱼ = zeros(NumDims,)
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

function getorthonormal(α, NumDims)
    # Initilize the matrix:
    A = Matrix(1.0 * I, NumDims, NumDims)
    A = reverse(A, dims=2)
    A[:, 1] = transpose(α)

    # Perform QR factorization:
    Q, _ = qr(A)
    Q = Matrix(Q)

    # Clean up the result:
    R = transpose(reverse(Q, dims=2))
    R = Matrix(R)

    return R
end