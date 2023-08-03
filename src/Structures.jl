#= Transformations =#
abstract type AbstractTransformation end

mutable struct NatafTransformation <: AbstractTransformation
    # Marginal distributions:
    X::Vector{<:Distribution}
    # Correlation matrix:
    ρˣ::Matrix{<:Real}
    # Distorted correlation matrix:
    ρᶻ::Matrix{Float64}
    # Lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix:
    L::LinearAlgebra.LowerTriangular{Float64,Matrix{Float64}}
    # Inverse of the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix:
    L⁻¹::LinearAlgebra.LowerTriangular{Float64,Matrix{Float64}}

    # Constructor:
    function NatafTransformation(X, ρˣ)
        ρᶻ, L, L⁻¹ = getdistortedcorrelation(X, ρˣ)
        new(X, ρˣ, ρᶻ, L, L⁻¹)
    end
end

mutable struct RosenblattTransformation <: AbstractTransformation

end

#= Sampling Techniques =#
abstract type AbstractSamplingTechnique end

# Inverse Transform Sampling:
struct ITS <: AbstractSamplingTechnique

end

# Latin Hypercube Sampling:
struct LHS <: AbstractSamplingTechnique

end

#= Reliability Analysis =#
abstract type AbstractReliabilityProblem end

abstract type AbstractReliabililyAnalysisMethod end
abstract type FORMSubmethod end
abstract type SORMSubmethod end

mutable struct ReliabilityProblem <: AbstractReliabilityProblem
    # Marginal random variables:
    X::Vector{<:Distribution}
    # Correlation matrix:
    ρˣ::Matrix{<:Real}
    # Limit state function:
    g::Function
end

# First-Order Reliability Method:
Base.@kwdef struct FORM <: AbstractReliabililyAnalysisMethod
    Submethod::FORMSubmethod = iHLRF()
end

struct MCFOSM <: FORMSubmethod # Mean-Centered First-Order Second-Moment method

end

struct MCFOSMCache
    β::Number
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
    β::Number
    PoF::Number
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
    β::Number
    PoF::Number
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
    β₁::Number
    PoF₁::Number
    β₂::Vector{Float64}
    PoF₂::Vector{Float64}
    H::Matrix{Float64}
    P::Matrix{Float64}
    A::Matrix{Float64}
    κ::Vector{Float64}
end

Base.@kwdef struct PF <: SORMSubmethod # Point-Fitting method

end

struct PFCache

end

