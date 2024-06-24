"""
    ReliabilityProblem <: AbstractReliabilityProblem

Type used to define reliability problems.

$(TYPEDFIELDS)
"""
mutable struct ReliabilityProblem <: AbstractReliabilityProblem
    "Random vector ``\\vec{X}``"
    X::AbstractVector{<:Distributions.UnivariateDistribution}
    "Correlation matrix ``\\rho^{X}``"
    ρˣ::AbstractMatrix{<:Real}
    "Limit state function ``g(\\vec{X})``"
    g::Function
end

include("MC.jl")
include("IS.jl")
include("FORM.jl")
include("SORM.jl")
include("SSM.jl")
include("DCM/DCM.jl")
include("DCM/GenzMalik.jl")