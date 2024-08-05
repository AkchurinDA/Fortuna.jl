"""
    FORM <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using First-Order Reliability Method (FORM).

$(TYPEDFIELDS)
"""
Base.@kwdef struct FORM <: AbstractReliabililyAnalysisMethod
    Submethod::FORMSubmethod = iHLRF()
end

"""
    MCFOSM <: FORMSubmethod

Type used to perform reliability analysis using Mean-Centered First-Order Second-Moment (MCFOSM) method.

$(TYPEDFIELDS)
"""
struct MCFOSM <: FORMSubmethod # Mean-Centered First-Order Second-Moment method

end

"""
    RF <: FORMSubmethod

Type used to perform reliability analysis using Rackwitz-Fiessler (RF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct RF <: FORMSubmethod # Rackwitz-Fiessler method
    "Maximum number of iterations"
    MaxNumIterations::Integer = 250
    "Convergance criterion ``\\epsilon``"
    ϵ::Real = 1E-6
end

"""
    HLRF <: FORMSubmethod

Type used to perform reliability analysis using Hasofer-Lind Rackwitz-Fiessler (HLRF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct HLRF <: FORMSubmethod # Hasofer-Lind Rackwitz-Fiessler method
    "Maximum number of iterations"
    MaxNumIterations::Integer = 250
    "Convergance criterion #1 ``\\epsilon_{1}``"
    ϵ₁::Real = 1E-6
    "Convergance criterion #1 ``\\epsilon_{2}``"
    ϵ₂::Real = 1E-6
    "Starting point ``x_{0}``"
    x₀::Union{Nothing, Vector{<:Real}} = nothing
end

"""
    iHLRF <: FORMSubmethod

Type used to perform reliability analysis using improved Hasofer-Lind Rackwitz-Fiessler (iHLRF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct iHLRF <: FORMSubmethod # Improved Hasofer-Lind Rackwitz-Fiessler method
    "Maximum number of iterations"
    MaxNumIterations::Integer = 250
    "Convergance criterion #1 ``\\epsilon_{1}``"
    ϵ₁::Real = 1E-6
    "Convergance criterion #1 ``\\epsilon_{2}``"
    ϵ₂::Real = 1E-6
    "Starting point ``x_{0}``"
    x₀::Union{Nothing, Vector{<:Real}} = nothing
    "c-coefficient applied at each iteration ``c_{0}``"
    c₀::Union{Nothing, Real} = nothing
end

"""
    MCFOSMCache

Type used to store results of reliability analysis performed using Mean-Centered First-Order Second-Moment (MCFOSM) method.

$(TYPEDFIELDS)
"""
struct MCFOSMCache
    "Reliability index ``\\beta``"
    β::Float64
end

"""
    RFCache

Type used to store results of reliability analysis performed using Rackwitz-Fiessler (RF) method.

$(TYPEDFIELDS)
"""
struct RFCache
    "Reliability index ``\\beta``"
    β::Float64
    "Design points in X-space at each iteration ``\\vec{x}_{i}^{*}``"
    x::Matrix{Float64}
    "Design points in U-space at each iteration ``\\vec{u}_{i}^{*}``"
    u::Matrix{Float64}
    "Means of equivalent normal marginals at each iteration ``\\vec{\\mu}_{i}``"
    μ::Matrix{Float64}
    "Standard deviations of equivalent normal marginals at each iteration ``\\vec{\\sigma}_{i}``"
    σ::Matrix{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla G(\\vec{u}_{i}^{*})``"
    ∇G::Matrix{Float64}
    "Normalized negative gradient of the limit state function at each iteration ``\\vec{\\alpha}_{i}``"
    α::Matrix{Float64}
    "Convergance status"
    Convergance::Bool
end

"""
    HLRFCache

Type used to store results of reliability analysis performed using Hasofer-Lind Rackwitz-Fiessler (HLRF) method.

$(TYPEDFIELDS)
"""
struct HLRFCache
    "Reliability index ``\\beta``"
    β::Float64
    "Probability of failure ``P_{f}``"
    PoF::Float64
    "Design points in X-space at each iteration ``\\vec{x}_{i}^{*}``"
    x::Matrix{Float64}
    "Design points in U-space at each iteration ``\\vec{u}_{i}^{*}``"
    u::Matrix{Float64}
    "Limit state function at each iteration ``G(\\vec{u}_{i}^{*})``"
    G::Vector{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla G(\\vec{u}_{i}^{*})``"
    ∇G::Matrix{Float64}
    "Normalized negative gradient of the limit state function at each iteration ``\\vec{\\alpha}_{i}``"
    α::Matrix{Float64}
    "Search direction at each iteration ``\\vec{d}_{i}``"
    d::Matrix{Float64}
    "Importance vector ``\\vec{\\gamma}``"
    γ::Vector{Float64}
    "Convergance status"
    Convergance::Bool
end

"""
    iHLRFCache

Type used to store results of reliability analysis performed using improved Hasofer-Lind Rackwitz-Fiessler (iHLRF) method.

$(TYPEDFIELDS)
"""
struct iHLRFCache
    "Reliability index ``\\beta``"
    β::Float64
    "Probability of failure ``P_{f}``"
    PoF::Float64
    "Design points in X-space at each iteration ``\\vec{x}_{i}^{*}``"
    x::Matrix{Float64}
    "Design points in U-space at each iteration ``\\vec{u}_{i}^{*}``"
    u::Matrix{Float64}
    "Limit state function at each iteration ``G(\\vec{u}_{i}^{*})``"
    G::Vector{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla G(\\vec{u}_{i}^{*})``"
    ∇G::Matrix{Float64}
    "Normalized negative gradient of the limit state function at each iteration ``\\vec{\\alpha}_{i}``"
    α::Matrix{Float64}
    "Search direction at each iteration ``\\vec{d}_{i}``"
    d::Matrix{Float64}
    "c-coefficient at each iteration ``c_{i}``"
    c::Vector{Float64}
    "Merit function at each iteration ``m_{i}``"
    m::Vector{Float64}
    "Step size at each iteration ``\\lambda_{i}``"
    λ::Vector{Float64}
    "Importance vector ``\\vec{\\gamma}``"
    γ::Vector{Float64}
    "Convergance status"
    Convergance::Bool
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::FORM; diff::Symbol = :automatic)

Function used to solve reliability problems using First-Order Reliability Method (FORM). \\
If `diff` is:
- `:automatic`, then the function will use automatic differentiation to compute gradients, jacobians, etc.
- `:numeric`, then the function will use numeric differentiation to compute gradients, jacobians, etc.
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::FORM; diff::Symbol = :automatic)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    # Extract the problem data:
    X  = Problem.X
    ρˣ = Problem.ρˣ
    g  = Problem.g

    if !isa(Submethod, MCFOSM) && !isa(Submethod, RF) && !isa(Submethod, HLRF) && !isa(Submethod, iHLRF)
        throw(ArgumentError("Invalid FORM submethod!"))
    elseif isa(Submethod, MCFOSM)
        # Compute the means of marginal distrbutions:
        Mˣ = Distributions.mean.(X)

        # Convert the correlation matrix into covariance matrix:
        σˣ = Distributions.std.(X)
        Dˣ = LinearAlgebra.diagm(σˣ)
        Σˣ = Dˣ * ρˣ * Dˣ

        # Compute gradient of the limit state function and evaluate it at the means of the marginal distributions:
        ∇g = if diff == :automatic
            try
                ForwardDiff.gradient(g, Mˣ)
            catch
                FiniteDiff.finite_difference_gradient(g, Mˣ)
            end
        elseif diff == :numeric
            FiniteDiff.finite_difference_gradient(g, Mˣ)
        end

        # Compute the reliability index:
        β = g(Mˣ) / sqrt(LinearAlgebra.transpose(∇g) * Σˣ * ∇g)

        # Return results:
        return MCFOSMCache(β)
    elseif isa(Submethod, RF)
        # Extract the analysis details:
        MaxNumIterations = Submethod.MaxNumIterations
        ϵ                = Submethod.ϵ

        # Error-catching:
        ρˣ == LinearAlgebra.I || throw(ArgumentError("RF method is only applicable to random vectors with uncorrelated marginals!"))

        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Preallocate:
        β  = Vector{Float64}(undef, MaxNumIterations)
        x  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        u  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        μ  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        σ  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        ∇G = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        α  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        Convergance = true

        # Initialize the design point in X-space:
        x[:, 1] = mean.(X)

        # Force the design point to lay on the failure boundary:
        function F(u, p)
            x′              = zeros(eltype(u), NumDimensions)
            x′[1:(end - 1)] = p[1:(end - 1)]
            x′[end]         = u
        
            return g(x′)
        end

        Problem   = NonlinearSolve.NonlinearProblem(F, mean(X[end]), x[:, 1])
        x[end, 1] = if diff == :automatic
            try
                Solution = NonlinearSolve.solve(Problem, nothing, abstol = 1E-9, reltol = 1E-9)
                Solution.u
            catch
                Solution = NonlinearSolve.solve(Problem, NonlinearSolve.FastShortcutNonlinearPolyalg(autodiff = NonlinearSolve.AutoFiniteDiff()), abstol = 1E-9, reltol = 1E-9)
                Solution.u
            end
        elseif diff == :numeric
            Solution = NonlinearSolve.solve(Problem, NonlinearSolve.FastShortcutNonlinearPolyalg(autodiff = NonlinearSolve.AutoFiniteDiff()), abstol = 1E-9, reltol = 1E-9)
            Solution.u
        end
        
        # Start iterating:
        for i in 1:MaxNumIterations
            # Compute the mean and standard deviation values of the equivalient normal marginals:
            for j in 1:NumDimensions
                σ[j, i] = Distributions.pdf(Distributions.Normal(), Distributions.quantile(Distributions.Normal(), Distributions.cdf(X[j], x[j, i]))) / Distributions.pdf(X[j], x[j, i])
                μ[j, i] = x[j, i] - σ[j, i] * Distributions.quantile(Distributions.Normal(), Distributions.cdf(X[j], x[j, i]))
            end

            # Compute the design point in U-space:
            u[:, i] = (x[:, i] - μ[:, i]) ./ σ[:, i]

            # Evaluate gradient of the limit state function at the design point in U-space:
            ∇G[:, i] = if diff == :automatic
                try
                    -σ[:, i] .* ForwardDiff.gradient(g, x[:, i])
                catch
                    -σ[:, i] .* FiniteDiff.finite_difference_gradient(g, x[:, i])
                end
            elseif diff == :numeric
                -σ[:, i] .* FiniteDiff.finite_difference_gradient(g, x[:, i])
            end

            # Compute the reliability index:
            β[i] = LinearAlgebra.dot(∇G[:, i], u[:, i]) / LinearAlgebra.norm(∇G[:, i])

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = ∇G[:, i] / LinearAlgebra.norm(∇G[:, i])

            # Check for convergance:
            if i != 1
                Criterion = abs(β[i] - β[i - 1])
                if Criterion < ϵ || i == MaxNumIterations
                    if i == MaxNumIterations
                        @warn """
                        RF method did not converge in the given maximum number of iterations (MaxNumIterations = $MaxNumIterations)!
                        Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criterion (ϵ)!
                        """

                        Convergance = false
                    end

                    # Clean up the results:
                    β  = β[i]
                    x  = x[:, 1:i]
                    u  = u[:, 1:i]
                    μ  = μ[:, 1:i]
                    σ  = σ[:, 1:i]
                    ∇G = ∇G[:, 1:i]
                    α  = α[:, 1:i]

                    # Return results:
                    return RFCache(β, x, u, μ, σ, ∇G, α, Convergance)

                    # Break out:
                    continue
                end
            end

            # Compute the new design point in U-space:
            u[:, i + 1] = β[i] * α[:, i]

            # Compute the new design point in X-space:
            x[:, i + 1] = μ[:, i] + σ[:, i] .* u[:, i + 1]

            # Force the design point to lay on the failure boundary:
            Problem       = NonlinearSolve.NonlinearProblem(F, x[end, i + 1], x[:, i + 1])
            x[end, i + 1] = if diff == :automatic
                try
                    Solution = NonlinearSolve.solve(Problem, nothing, abstol = 1E-9, reltol = 1E-9)
                    Solution.u
                catch
                    Solution = NonlinearSolve.solve(Problem, NonlinearSolve.FastShortcutNonlinearPolyalg(autodiff = NonlinearSolve.AutoFiniteDiff()), abstol = 1E-9, reltol = 1E-9)
                    Solution.u
                end
            elseif diff == :numeric
                Solution = NonlinearSolve.solve(Problem, NonlinearSolve.FastShortcutNonlinearPolyalg(autodiff = NonlinearSolve.AutoFiniteDiff()), abstol = 1E-9, reltol = 1E-9)
                Solution.u
            end
        end
    elseif isa(Submethod, HLRF)
        # Extract the analysis details:
        MaxNumIterations = Submethod.MaxNumIterations
        ϵ₁ = Submethod.ϵ₁
        ϵ₂ = Submethod.ϵ₂
        x₀ = Submethod.x₀

        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Preallocate:
        x  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        u  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        G  = Vector{Float64}(undef, MaxNumIterations)
        ∇G = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        α  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        d  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        Convergance = true

        # Perform the Nataf Transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Initialize the design point in X-space:
        x[:, 1] = isnothing(x₀) ? mean.(X) : x₀

        # Compute the initial design point in U-space:
        u[:, 1] = transformsamples(NatafObject, x[:, 1], :X2U)

        # Evaluate the limit state function at the initial design point:
        G₀ = g(x[:, 1])

        # Set the step size to unity:
        λ = 1

        # Start iterating:
        for i in 1:MaxNumIterations
            # Compute the Jacobian of the transformation of the design point from X- to U-space:
            Jₓᵤ = getjacobian(NatafObject, x[:, i], :X2U)

            # Evaluate the limit state function at the design point in X-space:
            G[i] = g(x[:, i])

            # Evaluate gradient of the limit state function at the design point in X-space:
            ∇g = if diff == :automatic
                try
                    LinearAlgebra.transpose(ForwardDiff.gradient(g, x[:, i]))
                catch
                    LinearAlgebra.transpose(FiniteDiff.finite_difference_gradient(g, x[:, i]))
                end
            elseif diff == :numeric
                LinearAlgebra.transpose(FiniteDiff.finite_difference_gradient(g, x[:, i]))
            end

            # Convert the evaluated gradient of the limit state function from X- to U-space:
            ∇G[:, i] = vec(∇g * Jₓᵤ)

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = -∇G[:, i] / LinearAlgebra.norm(∇G[:, i])

            # Compute the search direction:
            d[:, i] = (G[i] / LinearAlgebra.norm(∇G[:, i]) + LinearAlgebra.dot(α[:, i], u[:, i])) * α[:, i] - u[:, i]

            # Check for convergance:
            Criterion₁ = abs(g(x[:, i]) / G₀) # Check if the limit state function is close to zero.
            Criterion₂ = LinearAlgebra.norm(u[:, i] - LinearAlgebra.dot(α[:, i], u[:, i]) * α[:, i]) # Check if the design point is on the failure boundary.
            if (Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂) || i == MaxNumIterations
                # Check for convergance:
                if i == MaxNumIterations  
                    @warn """
                    HL-RF method did not converge in the given maximum number of iterations (MaxNumIterations = $MaxNumIterations)!
                    Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criteria (ϵ₁, ϵ₂)!
                    """

                    Convergance = false
                end

                # Compute the reliability index:
                β = LinearAlgebra.dot(α[:, i], u[:, i])

                # Compute the probability of failure:
                PoF = Distributions.cdf(Distributions.Normal(), -β)

                # Compute the importance vector:
                L⁻¹ = NatafObject.L⁻¹
                γ   = vec((LinearAlgebra.transpose(α[:, i]) * L⁻¹) / LinearAlgebra.norm(LinearAlgebra.transpose(α[:, i]) * L⁻¹))

                # Clean up the results:
                x  = x[:, 1:i]
                u  = u[:, 1:i]
                G  = G[1:i]
                ∇G = ∇G[:, 1:i]
                α  = α[:, 1:i]
                d  = d[:, 1:i]

                # Return results:
                return HLRFCache(β, PoF, x, u, G, ∇G, α, d, γ, Convergance)

                # Break out:
                continue
            end

            # Compute the new design point in U-space:
            u[:, i + 1] = u[:, i] + λ * d[:, i]

            # Compute the new design point in X-space:
            x[:, i + 1] = transformsamples(NatafObject, u[:, i + 1], :U2X)
        end
    elseif isa(Submethod, iHLRF)
        # Extract the analysis details:
        MaxNumIterations = Submethod.MaxNumIterations
        ϵ₁ = Submethod.ϵ₁
        ϵ₂ = Submethod.ϵ₂
        x₀ = Submethod.x₀
        c₀ = Submethod.c₀

        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Preallocate:
        x  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        u  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        G  = Vector{Float64}(undef, MaxNumIterations)
        ∇G = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        α  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        d  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        c  = Vector{Float64}(undef, MaxNumIterations - 1)
        m  = Vector{Float64}(undef, MaxNumIterations - 1)
        λ  = Vector{Float64}(undef, MaxNumIterations - 1)
        Convergance = true

        # Perform the Nataf Transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Initialize the design point in X-space:
        x[:, 1] = isnothing(x₀) ? mean.(X) : x₀

        # Compute the initial design point in U-space:
        u[:, 1] = transformsamples(NatafObject, x[:, 1], :X2U)

        # Evaluate the limit state function at the initial design point:
        G₀ = g(x[:, 1])

        # Start iterating:
        for i in 1:MaxNumIterations
            # Compute the Jacobian of the transformation of the design point from X- to U-space:
            Jₓᵤ = getjacobian(NatafObject, x[:, i], :X2U)

            # Evaluate the limit state function at the design point in X-space:
            G[i] = g(x[:, i])

            # Evaluate gradient of the limit state function at the design point in X-space:
            ∇g = if diff == :automatic
                try
                    LinearAlgebra.transpose(ForwardDiff.gradient(g, x[:, i]))
                catch
                    @warn "Automatic differentiation has failed. Switching to numeric differentiation."
                    LinearAlgebra.transpose(FiniteDiff.finite_difference_gradient(g, x[:, i]))
                end
            elseif diff == :numeric
                LinearAlgebra.transpose(FiniteDiff.finite_difference_gradient(g, x[:, i]))
            end

            # Convert the evaluated gradient of the limit state function from X- to U-space:
            ∇G[:, i] = vec(∇g * Jₓᵤ)

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = -∇G[:, i] / LinearAlgebra.norm(∇G[:, i])

            # Compute the search direction:
            d[:, i] = (G[i] / LinearAlgebra.norm(∇G[:, i]) + LinearAlgebra.dot(α[:, i], u[:, i])) * α[:, i] - u[:, i]

            # Check for convergance:
            Criterion₁ = abs(g(x[:, i]) / G₀) # Check if the limit state function is close to zero.
            Criterion₂ = LinearAlgebra.norm(u[:, i] - LinearAlgebra.dot(α[:, i], u[:, i]) * α[:, i]) # Check if the design point is on the failure boundary.
            if (Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂) || i == MaxNumIterations
                # Check for convergance:
                if i == MaxNumIterations  
                    @warn """
                    iHL-RF method did not converge in the given maximum number of iterations (MaxNumIterations = $MaxNumIterations)!
                    Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criteria (ϵ₁, ϵ₂)!
                    """

                    Convergance = false
                end

                # Compute the reliability index:
                β = LinearAlgebra.dot(α[:, i], u[:, i])

                # Compute the probability of failure:
                PoF = Distributions.cdf(Distributions.Normal(), -β)

                # Compute the importance vector:
                L⁻¹ = NatafObject.L⁻¹
                γ   = vec((LinearAlgebra.transpose(α[:, i]) * L⁻¹) / LinearAlgebra.norm(LinearAlgebra.transpose(α[:, i]) * L⁻¹))

                # Clean up the results:
                x  = x[:, 1:i]
                u  = u[:, 1:i]
                G  = G[1:i]
                ∇G = ∇G[:, 1:i]
                α  = α[:, 1:i]
                d  = d[:, 1:i]
                c  = c[1:(i - 1)]
                m  = m[1:(i - 1)]
                λ  = λ[1:(i - 1)]

                # Return results:
                return iHLRFCache(β, PoF, x, u, G, ∇G, α, d, c, m, λ, γ, Convergance)

                # Break out:
                continue
            end

            # Compute the c-coefficient:
            c[i] = isnothing(c₀) ? 2 * LinearAlgebra.norm(u[:, i]) / LinearAlgebra.norm(∇G[:, i]) + 10 : c₀

            # Compute the merit function at the current design point:
            m[i] = 0.5 * LinearAlgebra.norm(u[:, i]) ^ 2 + c[i] * abs(G[i])
            
            # Find a step size that satisfies m(uᵢ + λᵢdᵢ) < m(uᵢ):
            λₜ = 1
            uₜ = u[:, i] + λₜ * d[:, i]
            xₜ = transformsamples(NatafObject, uₜ, :U2X)
            Gₜ = g(xₜ)
            mₜ = 0.5 * LinearAlgebra.norm(uₜ) ^ 2 + c[i] * abs(Gₜ)
            ReductionCounter = 1
            while !(mₜ ≤ m[i])
                if ReductionCounter == 30
                    break
                end

                # Update the step size:
                λₜ = λₜ / 2

                # Recalculate the merit function:
                uₜ = u[:, i] + λₜ * d[:, i]
                xₜ = transformsamples(NatafObject, uₜ, :U2X)
                Gₜ = g(xₜ)
                mₜ = 0.5 * LinearAlgebra.norm(uₜ) ^ 2 + c[i] * abs(Gₜ)
                ReductionCounter = ReductionCounter + 1
            end

            # Update the step size:
            λ[i] = λₜ

            # Compute the new design point in U-space:
            u[:, i + 1] = u[:, i] + λ[i] * d[:, i]

            # Compute the new design point in X-space:
            x[:, i + 1] = transformsamples(NatafObject, u[:, i + 1], :U2X)
        end
    end
end