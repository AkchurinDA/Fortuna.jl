#=
Reliability analysis:
=#
mutable struct ReliabilityProblem
    # Marginal distributions:
    X::Vector{<:Distribution}
    # Correlation matrix:
    ρˣ::Matrix{<:Real}
    # Limit state function:
    g::Function
end

abstract type ReliabililyAnalysisMethod end

# Mean-Centered First-Order Second-Moment Method:
struct MCFOSM <: ReliabililyAnalysisMethod

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

struct FORM <: ReliabililyAnalysisMethod
    Submethod::FORMSubmethod
end

#=
Transformations:
=#
abstract type Transformation end

mutable struct NatafTransformation <: Transformation
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

mutable struct RosenblattTransformation <: Transformation

end