#=
Methods of reliability analysis:
=#
mutable struct ReliabilityProblem
    # Marginal distributions:
    X::Vector{<:Distribution}
    # Correlation matrix:
    ρˣ::Matrix{<:Real}
    # Limit state function:
    g::Function
end

#=
Transformations:
=#
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

mutable struct RosenblattTransformations <: AbstractTransformation

end