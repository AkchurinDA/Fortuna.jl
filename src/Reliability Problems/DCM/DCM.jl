"""
    DCM <: AbstractReliabililyAnalysisMethod

Type used to perform reliability analysis using the Divide-and-Conquer Method (DCM).

$(TYPEDFIELDS)
"""
Base.@kwdef struct DCM <: AbstractReliabililyAnalysisMethod
    "Number of initial divisions for each dimension"
    NumInitialDivisions::AbstractVector{<:Integer}
    "Maximum number of integration box subdivisions"
    MaxNumSubdivisions::Integer = 1E4
    "Relative tolerance for the error of the integration"
    ϵ::Real = 0.1
end

"""
    IntegrationBox

Type used to store the results of the integration.

$(TYPEDFIELDS)
"""
struct IntegrationBox
    a ::Vector{Float64}
    b ::Vector{Float64}
    IH::Float64
    IL::Float64
    E ::Float64
    K ::Integer
end

struct DCMCache
    I::Float64
    E::Float64
    # IntegrationBoxes::Vector{IntegrationBox}
    Convergance::Bool
end

function solve(Problem::ReliabilityProblem, AnalysisMethod::DCM; 
    RL::Union{Nothing, Real} = nothing, RU::Union{Nothing, Real} = nothing)
    # Extract analysis details:
    NumInitialDivisions = AnalysisMethod.NumInitialDivisions
    MaxNumSubdivisions  = AnalysisMethod.MaxNumSubdivisions
    ϵ                   = AnalysisMethod.ϵ

    # Extract the problem data:
    X  = Problem.X
    ρˣ = Problem.ρˣ
    g  = Problem.g

    # Perform the Nataf Transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Compute the number of dimensions: 
    NumDimensions = length(X)

    # Compute the bounds of the initial set of the integration boxes:
    # NOTE: There is definitely a huge potential for optimization here.
    RL = isnothing(RL) ?  0 : RL
    RU = isnothing(RU) ? 25 : RU
    TL = RL / (1 + RL)
    TU = RU / (1 + RU)
    if NumDimensions == 2
        LowerBounds = vcat(TL,     0 * ones(StaticArrays.SVector{NumDimensions - 1}))
        UpperBounds = vcat(TU, 2 * π * ones(StaticArrays.SVector{NumDimensions - 1}))
    else
        LowerBounds = vcat(TL, 0 * ones(StaticArrays.SVector{NumDimensions - 2}),     0)
        UpperBounds = vcat(TU, π * ones(StaticArrays.SVector{NumDimensions - 2}), 2 * π)
    end
    Δ = (UpperBounds - LowerBounds) ./ NumInitialDivisions
    a = [(LowerBounds[i]       ) .+ Δ[i] * range(0, NumInitialDivisions[i] - 1) for i in 1:NumDimensions]
    b = [(LowerBounds[i] + Δ[i]) .+ Δ[i] * range(0, NumInitialDivisions[i] - 1) for i in 1:NumDimensions]
    a = vec(collect.(Base.Iterators.product(a...)))
    b = vec(collect.(Base.Iterators.product(b...)))
    for i in 1:prod(NumInitialDivisions)
        a[i][1] = a[i][1] / (1  - a[i][1])
        b[i][1] = b[i][1] / (1  - b[i][1])
    end

    # Define the cubature rule to use:
    Rule = GenzMalikCubatureRule(Val{NumDimensions}())

    # Generate the initial set of the integration boxes:
    IntegrationBoxes = Vector{IntegrationBox}(undef, prod(NumInitialDivisions))
    ITotal = 0
    ETotal = 0
    for i in 1:prod(NumInitialDivisions)
        IntegrationBoxes[i] = IntegrationBox(a[i], b[i],  Rule(ξ -> integrand(ξ, g, NatafObject), a[i], b[i])...)

        ITotal += IntegrationBoxes[i].IH
        ETotal += IntegrationBoxes[i].E
    end

    # Check the convergance criterion:
    ETotal / LinearAlgebra.abs(ITotal) ≤ ϵ && return ITotal, ETotal

    for i in 1:MaxNumSubdivisions
        display(i)
        # Find the integration box with the largest error:
        sort!(IntegrationBoxes, by = x -> x.E)
        Box = pop!(IntegrationBoxes)

        # Subdivide the integration box along the dimension with the largest error:
        Δ = (Box.b[Box.K] - Box.a[Box.K]) / 2
        aNew = copy(Box.a)
        bNew = copy(Box.b)
        aNew[Box.K] += Δ
        bNew[Box.K] -= Δ

        Box1 = IntegrationBox(aNew, Box.b,  Rule(ξ -> integrand(ξ, g, NatafObject), aNew, Box.b)...)
        Box2 = IntegrationBox(Box.a, bNew,  Rule(ξ -> integrand(ξ, g, NatafObject), Box.a, bNew)...)

        push!(IntegrationBoxes, Box1)
        push!(IntegrationBoxes, Box2)

        ITotal += Box1.IH + Box2.IH - Box.IH
        ETotal += Box1.E  + Box2.E  - Box.E

        # Check the convergance criterion:
        if ETotal / LinearAlgebra.abs(ITotal) ≤ ϵ 
            println("Tolerance criterion is satisfied!")
            # return DCMCache(ITotal, ETotal, IntegrationBoxes, true)
            return DCMCache(ITotal, ETotal, true)
        end

        if i == MaxNumSubdivisions
            println("Maximum number of subdivisions is reached!")
            # return DCMCache(ITotal, ETotal, IntegrationBoxes, false)
            return DCMCache(ITotal, ETotal, false)
        end
    end
end

function integrand(ξ::AbstractVector{<:Real}, g::Function, NatafObject::NatafTransformation)
    # Compute the number of dimensions:
    NumDimensions = length(ξ)

    # Convert the coordinates from the spherical to Cartesian coordinate system:
    u = sphericaltocartesian(ξ)

    # Convert the coordinates from the standard normal space to the original space:
    x = transformsamples(NatafObject, u, :U2X)

    # Compute the indicator function:
    I = g(x) ≤ 0 ? 1 : 0

    # Compute the constrant term:
    if I == 1
        # Extract the spherical coordinates:
        r  = ξ[1]
        φ  = ξ[2:(end - 1)]

        # Compute the remaining terms of the integrand:
        C  = 1 / (sqrt(2 * π) ^ NumDimensions)
        Ψ₁ = (r ^ (NumDimensions - 1)) * exp(-r ^ 2 / 2)
        Ψ₂ = NumDimensions == 2 ? 1 : prod([sin(φ[i]) ^ (NumDimensions - 1 - i) for i in eachindex(φ)])

        # Return the results:
        return C * Ψ₁ * Ψ₂
    else
        # Return the results:
        return 0
    end 
end

function sphericaltocartesian(s::AbstractVector{<:Real})
    # Compute the number of dimensions:
    NumDimensions = length(s)

    x = Vector{Float64}(undef, NumDimensions)
    if NumDimensions == 2
        r = s[1]
        φ = s[2]

        x[1] = r * cos(φ)
        x[2] = r * sin(φ)
    else
        r = s[1]
        φ = s[2:end]

        for i in 1:NumDimensions
            if i == 1
                x[i] = r * cos(φ[i])
            elseif i == NumDimensions
                x[i] = r * prod(x -> sin(x), φ)
            else
                x[i] = r * prod(x -> sin(x), φ[1:(i - 1)]) * cos(φ[i])
            end
        end
    end

    return x    
end