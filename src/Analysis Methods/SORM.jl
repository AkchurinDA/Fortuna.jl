# Second-Order Reliability Method:
function analyze(Problem::ReliabilityProblem, AnalysisMethod::SORM)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    # Determine the design point using FORM:
    β₁, PoF₁, x, u = analyze(Problem, FORM(iHLRF()))
    x = x[:, end]
    u = u[:, end]

    if !isa(Submethod, CF) && !isa(Submethod, PF)
        error("Invalid SORM submethod.")
    elseif isa(Submethod, CF)
        # Extract the analysis details:
        ϵ = Submethod.ϵ

        # Extract data:
        g = Problem.g
        X = Problem.X
        ρˣ = Problem.ρˣ

        # Compute the number of marginal distributions:
        NumDims = length(X)

        # Perform Nataf transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Compute the Jacobian of the transformation of the design point from X- to U-space:
        Jₓᵤ = getjacobian(NatafObject, x, "X2U")

        # Evaluate gradient of the limit state function at the design point in X-space:
        ∇g = transpose(gradient(g, x))

        # Convert the evaluated gradient of the limit state function from X- to U-space:
        ∇G = ∇g * Jₓᵤ

        # Compute the normalized negative gradient vector at the design point in U-space:
        α = -∇G / norm(∇G)

        # Compute the Hessian at the design point in U-space:
        Hᵘ = gethessian(NatafObject, g, u, ϵ, NumDims)

        # Compute the orthonomal matrix:
        P = getorthonormal(α, NumDims)

        # Evaluate the principal curvatures:
        A = P * Hᵘ * transpose(P) / norm(∇G)
        κ = eigen(A[1:end-1, 1:end-1]).values

        # Compute the probabilities of failure:
        PoF₂ = Vector{Float64}(undef, 2)
        if all(x -> β₁ * x > -1, κ)
            # Hohenbichler and Rackwitz:
            ψ = pdf(Normal(0, 1), β₁) / cdf(Normal(0, 1), -β₁)
            PoF₂[1] = cdf(Normal(0, 1), -β₁) * prod(1 ./ sqrt.(1 .+ ψ .* κ))

            # Breaitung:
            PoF₂[2] = cdf(Normal(0, 1), -β₁) * prod(1 ./ sqrt.(1 .+ β₁ .* κ))

            # Generalized reliability index:
            β₂ = -quantile.(Normal(0, 1), PoF₂)
        else
            error("Approximations of the probability of failure from SORM are not valid since the principal curvatures are not located at the design point.")
        end

        # Return results:
        return β₁, β₂, PoF₁, PoF₂, κ
    elseif isa(Submethod, PF)

    end
end

function getorthonormal(α, NumDims)
    # Initilize the matrix:
    A = Matrix(1.0 * I, NumDims, NumDims)
    A = reverse(A, dims=2)
    A[:, 1] = transpose(α)

    # Perform QR factorization:
    Q, _ = qr(A)

    # Clean up the result:
    P = transpose(reverse(Q, dims=2))

    return P
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