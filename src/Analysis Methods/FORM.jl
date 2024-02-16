# First-Order Reliability Method:
"""
    analyze(Problem::ReliabilityProblem, AnalysisMethod::FORM)

The function solves the provided reliability problem using any submethod that falls under a broader category of First-Order Reliability Methods (FORM).
"""
function analyze(Problem::ReliabilityProblem, AnalysisMethod::FORM)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    # Extract the problem data:
    g = Problem.g
    X = Problem.X
    ρˣ = Problem.ρˣ

    if !isa(Submethod, MCFOSM) && !isa(Submethod, HL) && !isa(Submethod, RF) && !isa(Submethod, HLRF) && !isa(Submethod, iHLRF)
        error("Invalid FORM submethod.")
    elseif isa(Submethod, MCFOSM)
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

        # Return results:
        return MCFOSMCache(β)
    elseif isa(Submethod, HL)
        # Not yet implemented
    elseif isa(Submethod, RF)
        # Not yet implemented
    elseif isa(Submethod, HLRF)
        # Extract the analysis details:
        MaxNumIterations = Submethod.MaxNumIterations
        ϵ₁ = Submethod.ϵ₁
        ϵ₂ = Submethod.ϵ₂

        # Compute the number of marginal distributions:
        NumDims = length(X)

        # Preallocate:
        x   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        u   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        G   = Vector{Float64}(undef, MaxNumIterations)
        ∇G  = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        α   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        d   = Matrix{Float64}(undef, NumDims, MaxNumIterations)

        # Perform the Nataf Transformation:
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
            G[i] = g(x[:, i])

            # Evaluate gradient of the limit state function at the design point in X-space:
            ∇g = transpose(gradient(g, x[:, i]))

            # Convert the evaluated gradient of the limit state function from X- to U-space:
            ∇G[:, i] = vec(∇g * Jₓᵤ)

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = -∇G[:, i] / norm(∇G[:, i])

            # Compute the search direction:
            d[:, i] = (G[i] / norm(∇G[:, i]) + dot(α[:, i], u[:, i])) * α[:, i] - u[:, i]

            # Compute the new design point in U-space:
            u[:, i+1] = u[:, i] + λ * d[:, i]

            # Check for convergance:
            Criterion₁ = abs(g(x[:, i]) / G₀)                               # Check if the value of the limit state function is close to zero.
            Criterion₂ = norm(u[:, i] - dot(α[:, i], u[:, i]) * α[:, i])    # Check if the design point is on the failure boundary.
            if Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂
                # Compute the reliability index:
                β = dot(α[:, i], u[:, i])

                # Compute the probability of failure:
                PoF = cdf(Normal(0, 1), -β)

                # Clean up the results:
                x = x[:, 1:i]
                u = u[:, 1:i]
                G = G[1:i]
                ∇G = ∇G[:, 1:i]
                α = α[:, 1:i]
                d = d[:, 1:i]

                # Return results:
                return HLRFCache(β, PoF, x, u, G, ∇G, α, d)

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
        x   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        u   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        G   = Vector{Float64}(undef, MaxNumIterations)
        ∇G  = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        α   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        d   = Matrix{Float64}(undef, NumDims, MaxNumIterations)
        c   = Vector{Float64}(undef, MaxNumIterations)
        m   = Vector{Float64}(undef, MaxNumIterations)
        λ   = Vector{Float64}(undef, MaxNumIterations)

        # Perform the Nataf Transformation:
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
            G[i] = g(x[:, i])

            # Evaluate gradient of the limit state function at the design point in X-space:
            ∇g = transpose(gradient(g, x[:, i]))

            # Convert the evaluated gradient of the limit state function from X- to U-space:
            ∇G[:, i] = vec(∇g * Jₓᵤ)

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = -∇G[:, i] / norm(∇G[:, i])

            # Compute the search direction:
            d[:, i] = (G[i] / norm(∇G[:, i]) + dot(α[:, i], u[:, i])) * α[:, i] - u[:, i]

            # Compute the c-coefficient:
            c[i] = ceil(norm(u[:, i]) / norm(∇G[:, i]))

            # Compute the merit function at the current design point:
            m[i] = 0.5 * norm(u[:, i])^2 + c[i] * abs(G[i])

            # Find a step size that satisfies m(uᵢ + λᵢdᵢ) < m(uᵢ):
            λₜ = 1
            uₜ = u[:, i] + λₜ * d[:, i]
            xₜ = transformsamples(NatafObject, uₜ, "U2X")
            Gₜ = g(xₜ)
            mₜ = 0.5 * norm(uₜ)^2 + c[i] * abs(Gₜ)
            while mₜ ≥ m[i]
                # Update the step size:
                λₜ = λₜ / 2

                # Update the merit function:
                m[i] = mₜ

                # Recalculate the merit function:
                uₜ = u[:, i] + λₜ * d[:, i]
                xₜ = transformsamples(NatafObject, uₜ, "U2X")
                Gₜ = g(xₜ)
                mₜ = 0.5 * norm(uₜ)^2 + c[i] * abs(Gₜ)
            end

            # Update the step size:
            λ[i] = λₜ

            # Compute the new design point in U-space:
            u[:, i+1] = u[:, i] + λ[i] * d[:, i]

            # Compute the new design point in X-space:
            x[:, i+1] = transformsamples(NatafObject, u[:, i+1], "U2X")

            # Check for convergance:
            Criterion₁ = abs(g(x[:, i]) / G₀)                               # Check if the limit state function is close to zero.
            Criterion₂ = norm(u[:, i] - dot(α[:, i], u[:, i]) * α[:, i])    # Check if the design point is on the failure boundary.
            if Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂
                # Compute the reliability index:
                β = dot(α[:, i], u[:, i])

                # Compute the probability of failure:
                PoF = cdf(Normal(0, 1), -β)

                # Clean up the results:
                x   = x[:, 1:i]
                u   = u[:, 1:i]
                G   = G[1:i]
                ∇G  = ∇G[:, 1:i]
                α   = α[:, 1:i]
                d   = d[:, 1:i]
                c   = c[1:i]
                m   = m[1:i]
                λ   = λ[1:i]

                # Return results:
                return iHLRFCache(β, PoF, x, u, G, ∇G, α, d, c, m, λ)

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