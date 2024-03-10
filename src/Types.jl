# --------------------------------------------------
# SAMPLING TECHNIQUES
# --------------------------------------------------
abstract type AbstractSamplingTechnique end

# Inverse Transform Sampling:
struct ITS <: AbstractSamplingTechnique end

# Latin Hypercube Sampling:
struct LHS <: AbstractSamplingTechnique end

# --------------------------------------------------
# ISOPROBABILISTIC TRANSFORMATIONS
# --------------------------------------------------
abstract type AbstractTransformation end

mutable struct NatafTransformation <: AbstractTransformation
    X   ::AbstractVector{<:Distributions.Sampleable}
    ρˣ  ::AbstractMatrix{<:Real}
    ρᶻ  ::AbstractMatrix{Float64}
    L   ::AbstractMatrix{Float64}
    L⁻¹ ::AbstractMatrix{Float64}

    function NatafTransformation(X::AbstractVector{<:Distributions.Sampleable}, ρˣ::AbstractMatrix{<:Real})
        # Compute the distorted correlation matrix:
        ρᶻ, L, L⁻¹ = getdistortedcorrelation(X, ρˣ)

        # Return the Nataf Transformation object with the computed distorted correlation matrix:
        return new(X, ρˣ, ρᶻ, L, L⁻¹)
    end
end
Base.broadcastable(x::NatafTransformation) = Ref(x)

mutable struct RosenblattTransformation <: AbstractTransformation

end
Base.broadcastable(x::RosenblattTransformation) = Ref(x)

# --------------------------------------------------
# RELIABILITY PROBLEMS
# --------------------------------------------------
abstract type AbstractReliabilityProblem end

mutable struct ReliabilityProblem <: AbstractReliabilityProblem
    X   ::AbstractVector{<:Distributions.Sampleable}
    ρˣ  ::AbstractMatrix{<:Real}
    g   ::Function
end

abstract type AbstractReliabililyAnalysisMethod end

# Monte Carlo Simulations:
Base.@kwdef struct MC <: AbstractReliabililyAnalysisMethod
    NumSamples          ::Integer = 10 ^ 6
    SamplingTechnique   ::AbstractSamplingTechnique = ITS()
end

struct MCCache
    Samples ::Matrix{Float64}
    PoF     ::Float64
end

Base.@kwdef struct IS <: AbstractReliabililyAnalysisMethod
    q                   ::Distributions.Sampleable
    NumSamples          ::Integer = 10 ^ 6
end

struct ISCache
    Samples ::Matrix{Float64}
    PoF     ::Float64
end

abstract type FORMSubmethod end

Base.@kwdef struct FORM <: AbstractReliabililyAnalysisMethod
    Submethod::FORMSubmethod = iHLRF()
end

struct MCFOSM <: FORMSubmethod # Mean-Centered First-Order Second-Moment method

end

Base.@kwdef struct HL <: FORMSubmethod # Hasofer-Lind method

end

Base.@kwdef struct RF <: FORMSubmethod # Rackwitz-Fiessler method

end

Base.@kwdef struct HLRF <: FORMSubmethod # Hasofer-Lind Rackwitz-Fiessler method
    MaxNumIterations    ::Integer = 250
    ϵ₁                  ::Real = 10^(-9)
    ϵ₂                  ::Real = 10^(-9)
    x₀                  ::Union{Nothing, Vector{<:Real}} = nothing
end

Base.@kwdef struct iHLRF <: FORMSubmethod # Improved Hasofer-Lind Rackwitz-Fiessler method
    MaxNumIterations    ::Integer = 250
    ϵ₁                  ::Real = 10^(-9)
    ϵ₂                  ::Real = 10^(-9)
    x₀                  ::Union{Nothing, Vector{<:Real}} = nothing
end

struct MCFOSMCache
    β::Float64
end

struct HLCache

end

struct RFCache

end

struct HLRFCache
    β   ::Float64
    PoF ::Float64
    x   ::Matrix{Float64}
    u   ::Matrix{Float64}
    G   ::Vector{Float64}
    ∇G  ::Matrix{Float64}
    α   ::Matrix{Float64}
    d   ::Matrix{Float64}
    γ   ::Vector{Float64}
end

struct iHLRFCache
    β   ::Float64
    PoF ::Float64
    x   ::Matrix{Float64}
    u   ::Matrix{Float64}
    G   ::Vector{Float64}
    ∇G  ::Matrix{Float64}
    α   ::Matrix{Float64}
    d   ::Matrix{Float64}
    c   ::Vector{Float64}
    m   ::Vector{Float64}
    λ   ::Vector{Float64}
    γ   ::Vector{Float64}
end

abstract type SORMSubmethod end

Base.@kwdef struct SORM <: AbstractReliabililyAnalysisMethod
    Submethod::SORMSubmethod = CF()
end

Base.@kwdef struct CF <: SORMSubmethod # Curve-Fitting method
    ϵ::Real = 1 / 1000
end

Base.@kwdef struct PF <: SORMSubmethod # Point-Fitting method

end

struct CFCache # Curve-Fitting method
    β₁      ::Float64
    PoF₁    ::Float64
    β₂      ::Vector{Float64}
    PoF₂    ::Vector{Float64}
    H       ::Matrix{Float64}
    R       ::Matrix{Float64}
    A       ::Matrix{Float64}
    κ       ::Vector{Float64}
end

struct PFCache # Point-Fitting method
    β₁              ::Float64
    PoF₁            ::Float64
    β₂              ::Vector{Float64}
    PoF₂            ::Vector{Float64}
    FittingPoints⁻  ::Matrix{Float64}
    FittingPoints⁺  ::Matrix{Float64}
    κ₁              ::Matrix{Float64}
    κ₂              ::Matrix{Float64}
end

Base.@kwdef struct SSM <: AbstractReliabililyAnalysisMethod
    P₀              ::Real = 0.10
    NumSamples      ::Integer = 10000
    MaxNumSubsets   ::Integer = 25
end

struct SSMCache
    XSamplesSubset  ::Vector{Matrix{Float64}}
    USamplesSubset  ::Vector{Matrix{Float64}}
    CSubset         ::Vector{Float64}
    PoFSubset       ::Vector{Float64}
    PoF             ::Float64
end

# --------------------------------------------------
# INVERSE RELIABILITY PROBLEMS
# --------------------------------------------------
mutable struct InverseReliabilityProblem <: AbstractReliabilityProblem
    X   ::AbstractVector{<:Distributions.Sampleable}
    ρˣ  ::AbstractMatrix{<:Real}
    g   ::Function
    β   ::Real
end

# --------------------------------------------------
# SENSITIVITY PROBLEMS
# --------------------------------------------------
mutable struct SensitivityProblem <: AbstractReliabilityProblem
    X   ::AbstractVector{<:Distributions.Sampleable}
    ρˣ  ::AbstractMatrix{<:Real}
    g   ::Function
    θ   ::AbstractVector{<:Real}
end

struct SensitivityProblemCache
    FORMSolution    ::Union{HLRFCache, iHLRFCache}
    ∇β              ::Vector{Float64}
    ∇PoF            ::Vector{Float64}
end