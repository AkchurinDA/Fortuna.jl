# First-Order Reliability Method:
function FORM(Problem::ReliabilityProblem; MaxNumIterations=100, ϵ₁=10^(-9), ϵ₂=10^(-9))
    # Extract data:
    g = Problem.g
    X = Problem.X
    ρˣ = Problem.ρˣ

    # Compute the number of marginal distributions:
    NumDims = length(X)

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Initialize the design point in X-space:
    x = Matrix{Float64}(undef, NumDims, MaxNumIterations)
    x[:, 1] = mean.(X)

    # Compute the initial design point in U-space:
    u = Matrix{Float64}(undef, NumDims, MaxNumIterations)
    u[:, 1] = transformsamples(NatafObject, x[:, 1], "X2U")

    # Evaluate the limit state function at the initial design point:
    G₀ = g(x[:, 1])

    # Start iterating:
    for i in 1:MaxNumIterations-1
        # Compute the design point in X-space:
        if i != 1
            x[:, i] = transformsamples(NatafObject, u[:, i], "U2X")
        end

        # Compute the Jacobian of the transformation of the design point from X- to U-space:
        Jₓᵤ = getjacobian(NatafObject, x[:, i], "X2U")

        # Evaluate the limit state function at the design point in X-space:
        G = g(x[:, i])

        # Evaluate gradient of the limit state function at the design point in X-space:
        ∇g = transpose(gradient(g, x[:, i]))

        # Convert the evaluated gradient of the limit state function from X- to U-space:
        ∇G = ∇g * Jₓᵤ

        # Compute the normalized negative gradient vector at the design point in U-space:
        α = -∇G / norm(∇G)

        # Compute the search direction:
        d = (G / norm(∇G) + α * u[:, i]) * transpose(α) - u[:, i]

        # Compute the c-coefficient:
        c = ceil(norm(u[:, i]) / norm(∇G))

        # Compute the merit function:
        m = 0.5 * norm(u[:, i])^2 + c * abs(G)

        # Assume the step size:
        λ = 1

        # Check if the assumed step size satisfies the Armijo rule, such that m(uᵢ + λᵢdᵢ) < m(uᵢ):
        uₜ = u[:, i] + λ * d
        xₜ = transformsamples(NatafObject, uₜ, "U2X")
        Gₜ = g(xₜ)
        mₜ = 0.5 * norm(uₜ)^2 + c * abs(Gₜ)
        while mₜ ≥ m
            # Update the step size and merit function:
            λ = λ / 2
            m = mₜ

            # Recalculate the merit function:
            uₜ = u[:, i] + λ * d
            xₜ = transformsamples(NatafObject, uₜ, "U2X")
            Gₜ = g(xₜ)
            mₜ = 0.5 * norm(uₜ)^2 + c * abs(Gₜ)
        end

        # Update the merit function:
        m = mₜ

        # Compute the new design point in U-space:
        u[:, i+1] = u[:, i] + λ * d

        # Compute the new design point in X-space:
        x[:, i+1] = transformsamples(NatafObject, u[:, i+1], "U2X")

        # Check for convergance:
        Criterion₁ = abs(g(x[:, i+1]) / G₀) # Check if the design point is on the failure boundary.
        Criterion₂ = norm(u[:, i+1] - α * u[:, i+1] * transpose(α)) # Check if the design point is on the failure boundary.
        if Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂
            # Compute the Jacobian of the transformation of the design point from X- to U-space:
            Jₓᵤ = getjacobian(NatafObject, x[:, i+1], "X2U")

            # Evaluate the limit state function at the design point in X-space:
            G = g(x[:, i])

            # Evaluate gradient of the limit state function at the design point in X-space:
            ∇g = transpose(gradient(g, x[:, i+1]))

            # Convert the evauated gradient of the limit state function from X- to U-space:
            ∇G = ∇g * Jₓᵤ

            # Compute the normalized negative gradient vector at the design point in U-space:
            α = -∇G / norm(∇G)

            # Compute the reliability index:
            β = α * u[:, i+1]

            # Clean up the results:
            x = x[:, 1:i+1]
            u = u[:, 1:i+1]

            # Return results:
            return β, x, u

            # Break out:
            continue
        else
            if i == MaxNumIterations - 1
                error("FORM did not converge. Try increasing the maximum number of iterations allowed.")
            end
        end
    end
end