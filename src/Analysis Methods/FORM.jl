# First-Order Reliability Method:
function analyze(Problem::ReliabilityProblem, AnalysisMethod::FORM)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    if !isa(Submethod, MCFOSM) && !isa(Submethod, HL) && !isa(Submethod, RF) && !isa(Submethod, HLRF) && !isa(Submethod, iHLRF)
        error("Invalid FORM submethod.")
    elseif isa(Submethod, MCFOSM)
        # Extract data:
        X = Problem.X
        ρˣ = Problem.ρˣ
        g = Problem.g

        # Compute the means of marginal distrbutions:
        Mˣ = mean.(X)

        # Convert the correlation matrix into covariance matrix:
        σˣ = std.(X)
        Dˣ = diagm(σˣ)
        Σˣ = Dˣ * ρˣ * Dˣ

        # Compute gradient of the limit state function and evaluate it at the means of the marginal distributions:
        ∇g = gradient(g, Mˣ)

        # Compute the reliability index:
        β = g(Mˣ) / sqrt(transpose(∇g) * Σˣ * ∇g)

        return β
    elseif isa(Submethod, HL)
        # Not yet implemented
    elseif isa(Submethod, RF)
        # Not yet implemented
    elseif isa(Submethod, HLRF)
        # Extract the analysis details:
        MaxNumIterations = Submethod.MaxNumIterations
        ϵ₁ = Submethod.ϵ₁
        ϵ₂ = Submethod.ϵ₂

        # Extract data:
        g = Problem.g
        X = Problem.X
        ρˣ = Problem.ρˣ

        # Compute the number of marginal distributions:
        NumDims = length(X)

        # Preallocate:
        x = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        u = Matrix{Float64}(undef, NumDims, MaxNumIterations)

        # Perform Nataf transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Initialize the design point in X-space:
        x[:, 1] = mean.(X)

        # Compute the initial design point in U-space:
        u[:, 1] = transformsamples(NatafObject, x[:, 1], "X2U")

        # Evaluate the limit state function at the initial design point:
        G₀ = g(x[:, 1])

        # Set the step size to unity:
        λ = 1

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

                # Compute the probability of failure:
                PoF = cdf(Normal(0, 1), -β)

                # Clean up the results:
                x = x[:, 1:i+1]
                u = u[:, 1:i+1]

                # Return results:
                return β, PoF, x, u

                # Break out:
                continue
            else
                if i == MaxNumIterations - 1
                    error("HL-RF did not converge. Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criterions (ϵ₁ and ϵ₂).")
                end
            end
        end
    elseif isa(Submethod, iHLRF)
        # Extract the analysis details:
        MaxNumIterations = Submethod.MaxNumIterations
        ϵ₁ = Submethod.ϵ₁
        ϵ₂ = Submethod.ϵ₂

        # Extract data:
        g = Problem.g
        X = Problem.X
        ρˣ = Problem.ρˣ

        # Compute the number of marginal distributions:
        NumDims = length(X)

        # Preallocate:
        x = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        u = Matrix{Float64}(undef, NumDims, MaxNumIterations)

        # Perform Nataf transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Initialize the design point in X-space:
        x[:, 1] = mean.(X)

        # Compute the initial design point in U-space:
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

            # Assume the step size:
            λ = 1

            # Compute the c-coefficient:
            c = ceil(norm(u[:, i]) / norm(∇G))

            # Compute the merit function:
            m = 0.5 * norm(u[:, i])^2 + c * abs(G)

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

                # Compute the probability of failure:
                PoF = cdf(Normal(0, 1), -β)

                # Clean up the results:
                x = x[:, 1:i+1]
                u = u[:, 1:i+1]

                # Return results:
                return β, PoF, x, u

                # Break out:
                continue
            else
                if i == MaxNumIterations - 1
                    error("iHL-RF did not converge. Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criterions (ϵ₁ and ϵ₂).")
                end
            end
        end
    end
end