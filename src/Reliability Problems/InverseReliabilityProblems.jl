mutable struct InverseReliabilityProblem <: AbstractReliabilityProblem
    X  ::AbstractVector{<:Distributions.Sampleable}
    ρˣ ::AbstractMatrix{<:Real}
    g  ::Function
    β  ::Real
end