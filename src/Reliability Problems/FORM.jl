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

Base.@kwdef struct HL <: FORMSubmethod # Hasofer-Lind method

end

Base.@kwdef struct RF <: FORMSubmethod # Rackwitz-Fiessler method

end

"""
    HLRF <: FORMSubmethod

Type used to perform reliability analysis using Hasofer-Lind Rackwitz-Fiessler (HLRF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct HLRF <: FORMSubmethod # Hasofer-Lind Rackwitz-Fiessler method
    "Maximum number of iterations"
    MaxNumIterations    ::Integer = 250
    "Convergance criterion #1 ``\\epsilon_{1}``"
    ϵ₁                  ::Real = 10^(-9)
    "Convergance criterion #1 ``\\epsilon_{2}``"
    ϵ₂                  ::Real = 10^(-9)
    "Starting point ``x_{0}``"
    x₀                  ::Union{Nothing, Vector{<:Real}} = nothing
end

"""
    iHLRF <: FORMSubmethod

Type used to perform reliability analysis using improved Hasofer-Lind Rackwitz-Fiessler (iHLRF) method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct iHLRF <: FORMSubmethod # Improved Hasofer-Lind Rackwitz-Fiessler method
    "Maximum number of iterations"
    MaxNumIterations    ::Integer = 250
    "Convergance criterion #1 ``\\epsilon_{1}``"
    ϵ₁                  ::Real = 10^(-9)
    "Convergance criterion #1 ``\\epsilon_{2}``"
    ϵ₂                  ::Real = 10^(-9)
    "Starting point ``x_{0}``"
    x₀                  ::Union{Nothing, Vector{<:Real}} = nothing
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

struct HLCache

end

struct RFCache

end

"""
    HLRFCache

Type used to store results of reliability analysis performed using Hasofer-Lind Rackwitz-Fiessler (HLRF) method.

$(TYPEDFIELDS)
"""
struct HLRFCache
    "Reliability index ``\\beta``"
    β   ::Float64
    "Probability of failure ``P_{f}``"
    PoF ::Float64
    "Design points in X-space at each iteration ``\\vec{x}_{i}^{*}``"
    x   ::Matrix{Float64}
    "Design points in U-space at each iteration ``\\vec{u}_{i}^{*}``"
    u   ::Matrix{Float64}
    "Limit state function at each iteration ``G(\\vec{u}_{i}^{*})``"
    G   ::Vector{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla G(\\vec{u}_{i}^{*})``"
    ∇G  ::Matrix{Float64}
    "Normalized negative gradient of the limit state function at each iteration ``\\vec{\\alpha}_{i}``"
    α   ::Matrix{Float64}
    "Search direction at each iteration ``\\vec{d}_{i}``"
    d   ::Matrix{Float64}
    "Importance vector ``\\vec{\\gamma}``"
    γ   ::Vector{Float64}
end

"""
    iHLRFCache

Type used to store results of reliability analysis performed using improved Hasofer-Lind Rackwitz-Fiessler (iHLRF) method.

$(TYPEDFIELDS)
"""
struct iHLRFCache
    "Reliability index ``\\beta``"
    β   ::Float64
    "Probability of failure ``P_{f}``"
    PoF ::Float64
    "Design points in X-space at each iteration ``\\vec{x}_{i}^{*}``"
    x   ::Matrix{Float64}
    "Design points in U-space at each iteration ``\\vec{u}_{i}^{*}``"
    u   ::Matrix{Float64}
    "Limit state function at each iteration ``G(\\vec{u}_{i}^{*})``"
    G   ::Vector{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla G(\\vec{u}_{i}^{*})``"
    ∇G  ::Matrix{Float64}
    "Normalized negative gradient of the limit state function at each iteration ``\\vec{\\alpha}_{i}``"
    α   ::Matrix{Float64}
    "Search direction at each iteration ``\\vec{d}_{i}``"
    d   ::Matrix{Float64}
    "c-coefficient at each iteration ``c_{i}``"
    c   ::Vector{Float64}
    "Merit function at each iteration ``m_{i}``"
    m   ::Vector{Float64}
    "Step size at each iteration ``\\lambda_{i}``"
    λ   ::Vector{Float64}
    "Importance vector ``\\vec{\\gamma}``"
    γ   ::Vector{Float64}
end

"""
    solve(Problem::ReliabilityProblem, AnalysisMethod::FORM)

Function used to solve reliability problems using First-Order Reliability Method (FORM).
"""
function solve(Problem::ReliabilityProblem, AnalysisMethod::FORM)
    # Extract the analysis method:
    Submethod = AnalysisMethod.Submethod

    # Extract the problem data:
    X   = Problem.X
    ρˣ  = Problem.ρˣ
    g   = Problem.g

    if !isa(Submethod, MCFOSM) && !isa(Submethod, HL) && !isa(Submethod, RF) && !isa(Submethod, HLRF) && !isa(Submethod, iHLRF)
        error("Invalid FORM submethod.")
    elseif isa(Submethod, MCFOSM)
        # Compute the means of marginal distrbutions:
        Mˣ = Distributions.mean.(X)

        # Convert the correlation matrix into covariance matrix:
        σˣ = Distributions.std.(X)
        Dˣ = LinearAlgebra.diagm(σˣ)
        Σˣ = Dˣ * ρˣ * Dˣ

        # Compute gradient of the limit state function and evaluate it at the means of the marginal distributions:
        ∇g = ForwardDiff.gradient(g, Mˣ)

        # Compute the reliability index:
        β = g(Mˣ) / sqrt(LinearAlgebra.transpose(∇g) * Σˣ * ∇g)

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
        x₀ = Submethod.x₀

        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Preallocate:
        x   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        u   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        G   = Vector{Float64}(undef, MaxNumIterations)
        ∇G  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        α   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        d   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)

        # Perform the Nataf Transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Initialize the design point in X-space:
        if isnothing(x₀) 
            x[:, 1] = mean.(X)
        else
            x[:, 1] = x₀
        end

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
            ∇g = LinearAlgebra.transpose(ForwardDiff.gradient(g, x[:, i]))

            # Convert the evaluated gradient of the limit state function from X- to U-space:
            ∇G[:, i] = vec(∇g * Jₓᵤ)

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = -∇G[:, i] / LinearAlgebra.norm(∇G[:, i])

            # Compute the search direction:
            d[:, i] = (G[i] / LinearAlgebra.norm(∇G[:, i]) + LinearAlgebra.dot(α[:, i], u[:, i])) * α[:, i] - u[:, i]

            # Compute the new design point in U-space:
            u[:, i + 1] = u[:, i] + λ * d[:, i]

            # Check for convergance:
            Criterion₁ = abs(g(x[:, i]) / G₀)                               # Check if the value of the limit state function is close to zero.
            Criterion₂ = LinearAlgebra.norm(u[:, i] - LinearAlgebra.dot(α[:, i], u[:, i]) * α[:, i])    # Check if the design point is on the failure boundary.
            if Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂
                # Compute the reliability index:
                β = LinearAlgebra.dot(α[:, i], u[:, i])

                # Compute the probability of failure:
                PoF = Distributions.cdf(Distributions.Normal(), -β)

                # Compute the importance vector:
                L⁻¹ = NatafObject.L⁻¹
                γ   = vec((LinearAlgebra.transpose(α[:, i]) * L⁻¹) / LinearAlgebra.norm(LinearAlgebra.transpose(α[:, i]) * L⁻¹))

                # Clean up the results:
                x   = x[:, 1:i]
                u   = u[:, 1:i]
                G   = G[1:i]
                ∇G  = ∇G[:, 1:i]
                α   = α[:, 1:i]
                d   = d[:, 1:i]

                # Return results:
                return HLRFCache(β, PoF, x, u, G, ∇G, α, d, γ)

                # Break out:
                continue
            else
                # Check for convergance:
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
        x₀ = Submethod.x₀

        # Compute number of dimensions: 
        NumDimensions = length(X)

        # Preallocate:
        x   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        u   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        G   = Vector{Float64}(undef, MaxNumIterations)
        ∇G  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        α   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        d   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
        c   = Vector{Float64}(undef, MaxNumIterations)
        m   = Vector{Float64}(undef, MaxNumIterations)
        λ   = Vector{Float64}(undef, MaxNumIterations)

        # Perform the Nataf Transformation:
        NatafObject = NatafTransformation(X, ρˣ)

        # Initialize the design point in X-space:
        if isnothing(x₀) 
            x[:, 1] = Distributions.mean.(X)
        else
            x[:, 1] = x₀
        end

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
            ∇g = LinearAlgebra.transpose(ForwardDiff.gradient(g, x[:, i]))

            # Convert the evaluated gradient of the limit state function from X- to U-space:
            ∇G[:, i] = vec(∇g * Jₓᵤ)

            # Compute the normalized negative gradient vector at the design point in U-space:
            α[:, i] = -∇G[:, i] / LinearAlgebra.norm(∇G[:, i])

            # Compute the search direction:
            d[:, i] = (G[i] / LinearAlgebra.norm(∇G[:, i]) + LinearAlgebra.dot(α[:, i], u[:, i])) * α[:, i] - u[:, i]

            # Compute the c-coefficient:
            c[i] = ceil(LinearAlgebra.norm(u[:, i]) / LinearAlgebra.norm(∇G[:, i]))

            # Compute the merit function at the current design point:
            m[i] = 0.5 * LinearAlgebra.norm(u[:, i])^2 + c[i] * abs(G[i])

            # Find a step size that satisfies m(uᵢ + λᵢdᵢ) < m(uᵢ):
            λₜ = 1
            uₜ = u[:, i] + λₜ * d[:, i]
            xₜ = transformsamples(NatafObject, uₜ, "U2X")
            Gₜ = g(xₜ)
            mₜ = 0.5 * LinearAlgebra.norm(uₜ)^2 + c[i] * abs(Gₜ)
            while mₜ ≥ m[i]
                # Update the step size:
                λₜ = λₜ / 2

                # Update the merit function:
                m[i] = mₜ

                # Recalculate the merit function:
                uₜ = u[:, i] + λₜ * d[:, i]
                xₜ = transformsamples(NatafObject, uₜ, "U2X")
                Gₜ = g(xₜ)
                mₜ = 0.5 * LinearAlgebra.norm(uₜ)^2 + c[i] * abs(Gₜ)
            end

            # Update the step size:
            λ[i] = λₜ

            # Compute the new design point in U-space:
            u[:, i + 1] = u[:, i] + λ[i] * d[:, i]

            # Compute the new design point in X-space:
            x[:, i + 1] = transformsamples(NatafObject, u[:, i + 1], "U2X")

            # Check for convergance:
            Criterion₁ = abs(g(x[:, i]) / G₀)                               # Check if the limit state function is close to zero.
            Criterion₂ = LinearAlgebra.norm(u[:, i] - LinearAlgebra.dot(α[:, i], u[:, i]) * α[:, i])    # Check if the design point is on the failure boundary.
            if Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂
                # Compute the reliability index:
                β = LinearAlgebra.dot(α[:, i], u[:, i])

                # Compute the probability of failure:
                PoF = Distributions.cdf(Distributions.Normal(), -β)

                # Compute the importance vector:
                L⁻¹ = NatafObject.L⁻¹
                γ   = vec((LinearAlgebra.transpose(α[:, i]) * L⁻¹) / LinearAlgebra.norm(LinearAlgebra.transpose(α[:, i]) * L⁻¹))

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
                return iHLRFCache(β, PoF, x, u, G, ∇G, α, d, c, m, λ, γ)

                # Break out:
                continue
            else
                # Check for convergance:
                if i == MaxNumIterations - 1
                    error("iHL-RF did not converge. Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criterions (ϵ₁ and ϵ₂).")
                end
            end
        end
    end
end