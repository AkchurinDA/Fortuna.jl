"""
    InverseReliabilityProblem <: AbstractReliabilityProblem

Type used to define inverse reliability problems.

$(TYPEDFIELDS)
"""
mutable struct InverseReliabilityProblem <: AbstractReliabilityProblem
    "Random vector ``\\vec{X}``"
    X::AbstractVector{<:Distributions.UnivariateDistribution}
    "Correlation matrix ``\\rho^{X}``"
    ρˣ ::AbstractMatrix{<:Real}
    "Limit state function ``g(\\vec{X}, \\theta)``"
    g  ::Function
    "Target reliability index ``\\beta_t``"
    β  ::Real
end

"""
    InverseReliabilityProblemCache

Type used to store results of inverse reliability analysis.

$(TYPEDFIELDS)
"""
struct InverseReliabilityProblemCache
    "Design points in X-space at each iteration ``\\vec{x}_{i}^{*}``"
    x::Matrix{Float64}
    "Design points in U-space at each iteration ``\\vec{u}_{i}^{*}``"
    u::Matrix{Float64}
    "Parameter of interest at each iteration ``\\theta_{i}``"
    θ::Vector{Float64}
    "Limit state function at each iteration ``G(\\vec{u}_{i}^{*}, \\theta_{i})``"
    G::Vector{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla_{\\vec{u}} G(\\vec{u}_{i}^{*}, \\theta_{i})``"
    ∇Gu::Matrix{Float64}
    "Gradient of the limit state function at each iteration ``\\nabla_{\\theta} G(\\vec{u}_{i}^{*}, \\theta_{i})``"
    ∇Gθ::Vector{Float64}
    "Normalized negative gradient of the limit state function at each iteration ``\\vec{\\alpha}_{i}``"
    α::Matrix{Float64}
    "Search direction for the design point in U-space at each iteration ``\\vec{d}_{u_{i}}``"
    du::Matrix{Float64}
    "Search direction for the parameter of interest at each iteration ``\\vec{d}_{u_{i}}``"
    dθ::Vector{Float64}
    "``c_{1}``-coefficients at each iteration ``c_{1_{i}}``"
    c₁::Vector{Float64}
    "``c_{2}``-coefficients at each iteration ``c_{2_{i}}``"
    c₂::Vector{Float64}
    "First merit function at each iteration ``m_{1_{i}}``"
    m₁::Vector{Float64}
    "Second merit function at each iteration ``m_{2_{i}}``"
    m₂::Vector{Float64}
    "Merit function at each iteration ``m_{i}``"
    m::Vector{Float64}
    "Step size at each iteration ``\\lambda_{i}``"
    λ::Vector{Float64}
end

"""
    solve(Problem::InverseReliabilityProblem, θ₀::Real; 
        x₀::Union{Nothing, Vector{<:Real}} = nothing, 
        MaxNumIterations = 250, ϵ₁ = 10E-6, ϵ₂ = 10E-6, ϵ₃ = 10E-3,
        backend = AutoForwardDiff())

Function used to solve inverse reliability problems.
"""
function solve(Problem::InverseReliabilityProblem, θ₀::Real; 
    MaxNumIterations = 250, ϵ₁ = 1E-6, ϵ₂ = 1E-6, ϵ₃ = 1E-6,
    x₀::Union{Nothing, Vector{<:Real}} = nothing, 
    c₀::Union{Nothing, Real} = nothing,
    backend = AutoForwardDiff())
    # Extract the problem data:
    X  = Problem.X
    ρˣ = Problem.ρˣ
    g  = Problem.g
    β  = Problem.β

    # Compute number of dimensions: 
    NumDimensions = length(X)

    # Preallocate:
    x   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
    u   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
    θ   = Vector{Float64}(undef, MaxNumIterations)
    G   = Vector{Float64}(undef, MaxNumIterations)
    ∇Gu = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
    ∇Gθ = Vector{Float64}(undef, MaxNumIterations)
    α   = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
    du  = Matrix{Float64}(undef, NumDimensions, MaxNumIterations)
    dθ  = Vector{Float64}(undef, MaxNumIterations)
    c₁  = Vector{Float64}(undef, MaxNumIterations)
    c₂  = Vector{Float64}(undef, MaxNumIterations)
    m₁  = Vector{Float64}(undef, MaxNumIterations)
    m₂  = Vector{Float64}(undef, MaxNumIterations)
    m   = Vector{Float64}(undef, MaxNumIterations)
    λ   = Vector{Float64}(undef, MaxNumIterations)

    # Perform the Nataf Transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Initialize the design point in X-space:
    x[:, 1] = isnothing(x₀) ? mean.(X) : x₀

    # Initialize the unknown parameter:
    θ[1] = θ₀

    # Compute the initial design point in U-space:
    u[:, 1] = transformsamples(NatafObject, x[:, 1], :X2U)

    # Evaluate the limit state function at the initial design point:
    G₀ = g(x[:, 1], θ[1])

    # Start iterating:
    for i in 1:(MaxNumIterations - 1)
        # Compute the design point in X-space:
        if i != 1
            x[:, i] = transformsamples(NatafObject, u[:, i], :U2X)
        end

        # Compute the Jacobian of the transformation of the design point from X- to U-space:
        Jₓᵤ = getjacobian(NatafObject, x[:, i], :X2U)

        # Evaluate the limit state function at the design point in X-space:
        G[i] = g(x[:, i], θ[i])

        # Evaluate gradients of the limit state function at the design point in X-space:
        ∇gx = try
            local ∇gx(x, θ) = LinearAlgebra.transpose(ForwardDiff.gradient(Unknown -> g(Unknown, θ), x))
            ∇gx(x[:, i], θ[i])
        catch
            local ∇gx(x, θ) = LinearAlgebra.transpose(FiniteDiff.finite_difference_gradient(Unknown -> g(Unknown, θ), x))
            ∇gx(x[:, i], θ[i])
        end

        ∇gθ = try
            local ∇gθ(x, θ) = ForwardDiff.derivative(Unknown -> g(x, Unknown), θ)
            ∇gθ(x[:, i], θ[i])
        catch
            local ∇gθ(x, θ) = FiniteDiff.finite_difference_derivative(Unknown -> g(x, Unknown), θ)
            ∇gθ(x[:, i], θ[i])
        end

        # Convert the evaluated gradients of the limit state function from X- to U-space:
        ∇Gu[:, i] = vec(∇gx * Jₓᵤ)
        ∇Gθ[i]    = ∇gθ

        # Compute the normalized negative gradient vector at the design point in U-space:
        α[:, i] = -∇Gu[:, i] / LinearAlgebra.norm(∇Gu[:, i])

        # Compute the c-coefficients:
        c₁[i] = isnothing(c₀) ? 2 * LinearAlgebra.norm(u[:, i]) / LinearAlgebra.norm(∇Gu[:, i]) + 10 : c₀
        c₂[i] = 1

        # Compute the merit functions at the current design point:
        m₁[i] = 0.5 * LinearAlgebra.norm(u[:, i]) ^ 2 + c₁[i] * abs(G[i])
        m₂[i] = 0.5 * c₂[i] * (LinearAlgebra.norm(u[:, i]) - β) ^ 2
        m[i]  = m₁[i] + m₂[i]

        # Compute the search directions:
        du[:, i] = β * α[:, i] - u[:, i]
        dθ[i]    = (LinearAlgebra.norm(∇Gu[:, i]) / ∇Gθ[i]) * (β - LinearAlgebra.dot(α[:, i], u[:, i]) - G[i] / LinearAlgebra.norm(∇Gu[:, i]))

        # Find a step size that satisfies m(uᵢ + λᵢdᵢ) < m(uᵢ):
        λₜ = 1
        uₜ = u[:, i] + λₜ * du[:, i]
        θₜ = θ[i] + λₜ * dθ[i]
        xₜ = transformsamples(NatafObject, uₜ, :U2X)
        Gₜ = g(xₜ, θₜ)
        m₁ₜ = 0.5 * LinearAlgebra.norm(uₜ) ^ 2 + c₁[i] * abs(Gₜ)
        m₂ₜ = 0.5 * c₂[i] * (LinearAlgebra.norm(u[:, i]) - β) ^ 2
        mₜ  = m₁ₜ + m₂ₜ
        while mₜ > m[i]
            # Update the step size:
            λₜ = λₜ / 2

            # Recalculate the merit function:
            uₜ = u[:, i] + λₜ * du[:, i]
            θₜ = θ[i] + λₜ * dθ[i]
            xₜ = transformsamples(NatafObject, uₜ, :U2X)
            Gₜ = g(xₜ, θₜ)
            m₁ₜ = 0.5 * LinearAlgebra.norm(uₜ) ^ 2 + c₁[i] * abs(Gₜ)
            m₂ₜ = 0.5 * c₂[i] * (LinearAlgebra.norm(uₜ) - β) ^ 2
            mₜ = m₁ₜ + m₂ₜ
        end

        # Update the step size:
        λ[i] = λₜ

        # Compute the new design point in U-space:
        u[:, i + 1] = u[:, i] + λ[i] * du[:, i]
        θ[i + 1] = θ[i] + λ[i] * dθ[i]

        # Compute the new design point in X-space:
        x[:, i + 1] = transformsamples(NatafObject, u[:, i + 1], :U2X)

        # Check for convergance:
        Criterion₁ = abs(g(x[:, i], θ[i]) / G₀) # Check if the limit state function is close to zero.
        Criterion₂ = LinearAlgebra.norm(u[:, i] - LinearAlgebra.dot(α[:, i], u[:, i]) * α[:, i]) # Check if the design point is on the failure boundary.
        Criterion₃ = LinearAlgebra.norm(u[:, i + 1] - u[:, i]) / LinearAlgebra.norm(u[:, i]) 
                   + abs(θ[i + 1] - θ[i]) / abs(θ[i]) 
                   + abs(LinearAlgebra.dot(α[:, i], u[:, i]) - β) / β # Check if the solution has converged.
        if Criterion₁ < ϵ₁ && Criterion₂ < ϵ₂ && Criterion₃ < ϵ₃ && i != MaxNumIterations
            # Clean up the results:
            x   = x[:, 1:i]
            u   = u[:, 1:i]
            θ   = θ[1:i]
            G   = G[1:i]
            ∇Gu = ∇Gu[:, 1:i]
            ∇Gθ = ∇Gθ[1:i]
            α   = α[:, 1:i]
            du  = du[:, 1:i]
            dθ  = dθ[1:i]
            c₁  = c₁[1:i]
            c₂  = c₂[1:i]
            m₁  = m₁[1:i]
            m₂  = m₂[1:i]
            m   = m[1:i]
            λ   = λ[1:i]
            
            # Return results:
            return InverseReliabilityProblemCache(x, u, θ, G, ∇Gu, ∇Gθ, α, du, dθ, c₁, c₂, m₁, m₂, m, λ)

            # Break out:
            continue
        else
            # Check for convergance:
            i == MaxNumIterations && error("The solution did not converge. Try increasing the maximum number of iterations (MaxNumIterations) or relaxing the convergance criterions (ϵ₁, ϵ₂, and ϵ₃).")
        end
    end
end