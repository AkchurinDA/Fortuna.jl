"""
    SensitivityProblemTypeI <: AbstractReliabilityProblem

Type used to define sensitivity problems of type I (sensitivities w.r.t. the parameters of the limit state function).

$(TYPEDFIELDS)
"""
mutable struct SensitivityProblemTypeI <: AbstractReliabilityProblem
    "Random vector ``\\vec{X}``"
    X::AbstractVector{<:Distributions.UnivariateDistribution}
    "Correlation matrix ``\\rho^{X}``"
    ρˣ::AbstractMatrix{<:Real}
    "Limit state function ``g(\\vec{X}, \\vec{\\Theta})``"
    g::Function
    "Parameters of limit state function ``\\vec{\\Theta}``"
    Θ::AbstractVector{<:Real}
end

"""
    SensitivityProblemTypeII <: AbstractReliabilityProblem

Type used to define sensitivity problems of type II (sensitivities w.r.t. the parameters of the random vector).

$(TYPEDFIELDS)
"""
mutable struct SensitivityProblemTypeII <: AbstractReliabilityProblem
    "Random vector ``\\vec{X}(\\vec{\\Theta})``"
    X::Function
    "Correlation matrix ``\\rho^{X}``"
    ρˣ::AbstractMatrix{<:Real}
    "Limit state function ``g(\\vec{X})``"
    g::Function
    "Parameters of limit state function ``\\vec{\\Theta}``"
    Θ::AbstractVector{<:Real}
end


"""
    SensitivityProblemCache

Type used to store results of sensitivity analysis for problems of type I (sensitivities w.r.t. the parameters of the limit state function).

$(TYPEDFIELDS)
"""
struct SensitivityProblemCache
    "Results of reliability analysis performed using First-Order Reliability Method (FORM)"
    FORMSolution::iHLRFCache
    "Sensivity vector of reliability index ``\\vec{\\nabla}_{\\vec{\\Theta}} \\beta``"
    ∇β::Vector{Float64}
    "Sensivity vector of probability of failure ``\\vec{\\nabla}_{\\vec{\\Theta}} P_{f}``"
    ∇PoF::Vector{Float64}
end

"""
    solve(Problem::SensitivityProblemTypeI; backend = AutoForwardDiff())

Function used to solve sensitivity problems of type I (sensitivities w.r.t. the parameters of the limit state function).
"""
function solve(Problem::SensitivityProblemTypeI; backend = AutoForwardDiff())
    # Extract the problem data:
    X  = Problem.X
    ρˣ = Problem.ρˣ
    g  = Problem.g
    Θ  = Problem.Θ

    # Define a reliability problem for the FORM analysis:
    g₁(x) = g(x, Θ)
    FORMProblem = ReliabilityProblem(X, ρˣ, g₁)

    # Solve the reliability problem using the FORM:
    FORMSolution = solve(FORMProblem, FORM(), backend = backend)
    x            = FORMSolution.x[:, end]
    u            = FORMSolution.u[:, end]
    β            = FORMSolution.β

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Define gradient functions of the limit state function in X- and U-spaces, and compute the sensitivities vector for the reliability index:
    ∇β = try
        local ∇g(x, θ) = gradient(Unknown -> g(x, Unknown), backend, θ)
        local ∇G(u, θ) = gradient(Unknown -> G(g, θ, NatafObject, Unknown), backend, u)
        ∇g(x, Θ) / LinearAlgebra.norm(∇G(u, Θ))
    catch
        local ∇g(x, θ) = gradient(Unknown -> g(x, Unknown), AutoFiniteDiff(), θ)
        local ∇G(u, θ) = gradient(Unknown -> G(g, θ, NatafObject, Unknown), AutoFiniteDiff(), u)
        ∇g(x, Θ) / LinearAlgebra.norm(∇G(u, Θ))
    end

    # Compute the sensitivities vector for the probability of failure:
    ∇PoF = -Distributions.pdf(Distributions.Normal(), β) * ∇β

    return SensitivityProblemCache(FORMSolution, ∇β, ∇PoF)
end

"""
    solve(Problem::SensitivityProblemTypeII; backend = AutoForwardDiff())

Function used to solve sensitivity problems of type II (sensitivities w.r.t. the parameters of the random vector).
"""
function solve(Problem::SensitivityProblemTypeII; backend = AutoForwardDiff())
    # Extract the problem data:
    X  = Problem.X
    ρˣ = Problem.ρˣ
    g  = Problem.g
    Θ  = Problem.Θ

    # Define a reliability problem for the FORM analysis:
    FORMProblem = ReliabilityProblem(X(Θ), ρˣ, g)

    # Solve the reliability problem using the FORM:
    FORMSolution = solve(FORMProblem, FORM(), backend = backend)
    x            = FORMSolution.x[:, end]
    α            = FORMSolution.α[:, end]
    β            = FORMSolution.β

    # Define the Jacobian of the transformation function w.r.t. the parameters of the random vector and compute the sensitivity vector for the reliability index:
    ∇β = try
        local ∇T(θ) = jacobian(Unknown -> transformsamples(NatafTransformation(X(Unknown), ρˣ), x, :X2U), backend, θ)
        vec(LinearAlgebra.transpose(α) * ∇T(Θ))
    catch
        local ∇T(θ) = jacobian(Unknown -> transformsamples(NatafTransformation(X(Unknown), ρˣ), x, :X2U), AutoFiniteDiff(), θ)
        vec(LinearAlgebra.transpose(α) * ∇T(Θ))
    end

    # Compute the sensitivity vector for the probability of failure:
    ∇PoF = -Distributions.pdf(Distributions.Normal(), β) * ∇β

    return SensitivityProblemCache(FORMSolution, ∇β, ∇PoF)
end

function G(g::Function, Θ::AbstractVector{<:Real}, NatafObject::NatafTransformation, USample::AbstractVector{<:Real})
    # Transform samples:
    XSample = transformsamples(NatafObject, USample, :U2X)

    # Evaluate the limit state function at the transform samples:
    GSample = g(XSample, Θ)

    # Return the result:
    return GSample
end