#= Transformations =#
"""
    AbstractTransformation

A custom abstract supertype used by `Fortuna.jl` to define various types of isoprobabilistic transformations.
"""
abstract type AbstractTransformation end

"""
    mutable struct NatafTransformation <: AbstractTransformation

A custom type used by `Fortuna.jl` to perform the Nataf Transformation.

$(TYPEDFIELDS)
"""
mutable struct NatafTransformation <: AbstractTransformation
    "Random vector with correlated non-normal marginal random variables"
    X::Vector{<:Distribution}
    "Correlation matrix of the random vector ``\\vec{X}``"
    ρˣ::Matrix{<:Real}
    "Distorted correlation matrix the random vector ``\\vec{Z}``"
    ρᶻ::Matrix{Float64}
    "Lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix ``\\rho^{Z}``"
    L::LinearAlgebra.LowerTriangular{Float64,Matrix{Float64}}
    "Inverse of the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix ``\\rho^{Z}``"
    L⁻¹::LinearAlgebra.LowerTriangular{Float64,Matrix{Float64}}

    @doc """
        NatafTransformation(X::Vector{<:Distribution}, ρˣ::Matrix{<:Real})

    An inner constructor used to create a `NatafTransformation` object.

    $(TYPEDFIELDS)
    """
    function NatafTransformation(X::Vector{<:Distribution}, ρˣ::Matrix{<:Real})
        ρᶻ, L, L⁻¹ = getdistortedcorrelation(X, ρˣ)
        new(X, ρˣ, ρᶻ, L, L⁻¹)
    end
end

"""
    mutable struct RosenblattTransformation <: AbstractTransformation

A custom type used by `Fortuna.jl` to perform the Rosenblatt Transformation.
"""
mutable struct RosenblattTransformation <: AbstractTransformation

end

#= Sampling Techniques =#
"""
    AbstractSamplingTechnique

A custom abstract supertype used by `Fortuna.jl` to define various types of sampling techniques.
"""
abstract type AbstractSamplingTechnique end

# Inverse Transform Sampling:
struct ITS <: AbstractSamplingTechnique

end

# Latin Hypercube Sampling:
struct LHS <: AbstractSamplingTechnique

end

#= Reliability Analysis =#
"""
    AbstractReliabilityProblem

A custom abstract supertype used by `Fortuna.jl` to define various types of reliability problems.
"""
abstract type AbstractReliabilityProblem end

"""
    AbstractReliabililyAnalysisMethod

A custom abstract supertype used by `Fortuna.jl` to define various types of reliability analysis methods.
"""
abstract type AbstractReliabililyAnalysisMethod end
abstract type FORMSubmethod end
abstract type SORMSubmethod end

"""
    mutable struct ReliabilityProblem <: AbstractReliabilityProblem

A custom type used by `Fortuna.jl` to define reliability problems.

$(TYPEDFIELDS)
"""
mutable struct ReliabilityProblem <: AbstractReliabilityProblem
    "Random vector with correlated non-normal marginal random variables"
    X::Vector{<:Distribution}
    "Correlation matrix"
    ρˣ::Matrix{<:Real}
    "Limit state function"
    g::Function
end

# First-Order Reliability Method:
Base.@kwdef struct FORM <: AbstractReliabililyAnalysisMethod
    Submethod::FORMSubmethod = iHLRF()
end

struct MCFOSM <: FORMSubmethod # Mean-Centered First-Order Second-Moment method

end

struct MCFOSMCache
    β::Float64
end

Base.@kwdef struct HL <: FORMSubmethod # Hasofer-Lind method

end

struct HLCache

end

Base.@kwdef struct RF <: FORMSubmethod # Rackwitz-Fiessler method

end

struct RFCache

end

Base.@kwdef struct HLRF <: FORMSubmethod # Hasofer-Lind Rackwitz-Fiessler method
    # Maximum number of iterations allowed:
    MaxNumIterations::Integer = 100
    # Criterion #1:
    ϵ₁::Number = 10^(-9)
    # Criterion #2:
    ϵ₂::Number = 10^(-9)
end

struct HLRFCache
    β::Float64
    PoF::Float64
    x::Matrix{Float64}
    u::Matrix{Float64}
    G::Vector{Float64}
    ∇G::Matrix{Float64}
    α::Matrix{Float64}
    d::Matrix{Float64}
end

Base.@kwdef struct iHLRF <: FORMSubmethod # Improved Hasofer-Lind Rackwitz-Fiessler method
    # Maximum number of iterations allowed:
    MaxNumIterations::Integer = 100
    # Criterion #1:
    ϵ₁::Number = 10^(-9)
    # Criterion #2:
    ϵ₂::Number = 10^(-9)
end

struct iHLRFCache
    β::Float64
    PoF::Float64
    x::Matrix{Float64}
    u::Matrix{Float64}
    G::Vector{Float64}
    ∇G::Matrix{Float64}
    α::Matrix{Float64}
    d::Matrix{Float64}
    c::Vector{Float64}
    m::Vector{Float64}
    λ::Vector{Float64}
end


# Second-Order Reliability Method:
Base.@kwdef struct SORM <: AbstractReliabililyAnalysisMethod
    Submethod::SORMSubmethod = CF()
end

Base.@kwdef struct CF <: SORMSubmethod # Curve-Fitting method
    ϵ::Number = 1 / 1000
end

struct CFCache
    β₁::Float64
    PoF₁::Float64
    β₂::Vector{Float64}
    PoF₂::Vector{Float64}
    H::Matrix{Float64}
    R::Matrix{Float64}
    A::Matrix{Float64}
    κ::Vector{Float64}
end

Base.@kwdef struct PF <: SORMSubmethod # Point-Fitting method

end

struct PFCache

end

# Subset Simulation Method:
"""
    Base.@kwdef struct SSM <: AbstractReliabililyAnalysisMethod

A custom type used by `Fortuna.jl` to perform the Subset Siumlation Method.

$(TYPEDFIELDS)
"""
Base.@kwdef struct SSM <: AbstractReliabililyAnalysisMethod
    "Target conditional probability"
    P₀::Float64 = 0.1
    "Number of samples to generate for each subset"
    NumSamples::Integer = 1000
    "Maximum number of subsets"
    MaxNumSubsets::Integer = 10
end

"""
    struct SSMCache

A custom type used by `Fortuna.jl` to store the results of the analysis performed using the Subset Siumlation Method.

$(TYPEDFIELDS)
"""
struct SSMCache
    "Samples generated within each subset in X-space"
    XSamplesSubset::Vector{Matrix{Float64}}
    "Samples generated within each subset in U-space"
    USamplesSubset::Vector{Matrix{Float64}}
    "Thresholds for each subset"
    CSubset::Vector{Float64}
    "Probabilities of failure for each subset"
    PoFSubset::Vector{Float64}
    "Probabilities of failure"
    PoF::Float64
end