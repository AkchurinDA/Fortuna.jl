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
mutable struct ReliabilityProblem
    # Marginal distributions:
    X::Vector{<:Distribution}
    # Correlation matrix:
    ρˣ::Matrix{<:Real}
    # Limit state function:
    g::Function
end

abstract type AbstractReliabililyAnalysisMethod end

# Mean-Centered First-Order Second-Moment Method:
struct MCFOSM <: AbstractReliabililyAnalysisMethod

end

# First-Order Reliability Method:
abstract type FORMSubmethod end

Base.@kwdef struct HLRF <: FORMSubmethod # Hasofer-Lind Rackwitz-Fiessler method

end

Base.@kwdef struct iHLRF <: FORMSubmethod # Improved Hasofer-Lind Rackwitz-Fiessler method
    MaxNumIterations::Integer = 100
    ϵ₁::Number = 10^(-9)
    ϵ₂::Number = 10^(-9)
end

struct FORM <: AbstractReliabililyAnalysisMethod
    Submethod::FORMSubmethod
end
