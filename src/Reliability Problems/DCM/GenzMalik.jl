struct GenzMalikCubatureRule{N}
    P ::NTuple{4, Vector{StaticArrays.SVector{N}}} # Integration points for the last 4 fully symmetric sum terms
    WH::NTuple{5}                                  # Weights 5-point Genz-Malik cubature rule
    WL::NTuple{4}                                  # Weights 4-point Genz-Malik cubature rule
end

# Cache the Genz-Malik cubature rule's integration locations and weights:
const GMCRCache = Dict{Integer, GenzMalikCubatureRule}()

# Define an outer constructor for the Genz-Malik cubature rule:
function GenzMalikCubatureRule(V::Val{N}) where {N}
    # Check if the Genz-Malik cubature rule has already been cached:
    haskey(GMCRCache, N) && return GMCRCache[N]::GenzMalikCubatureRule{N}

    # Error-catching:
    N < 2 && throw(ArgumentError("Genz-Malik cubature rule works only for multi-dimensional integrals (N ≥ 2)!"))

    # Compute the integration points:
    λ₂ = sqrt(9 / 70)
    λ₃ = sqrt(9 / 10)
    λ₄ = sqrt(9 / 10)
    λ₅ = sqrt(9 / 19)

    # Compute all permutations of the integration points:
    P₂ = combinations(λ₂, V, 1) # In the paper, these should also include sign combinations. However, for the ease of computing integration errors along each dimensions these do not include sign combination.
    P₃ = combinations(λ₃, V, 1) # In the paper, these should also include sign combinations. However, for the ease of computing integration errors along each dimensions these do not include sign combination.
    P₄ = signedcombinations(λ₄, V, 2)
    P₅ = signedcombinations(λ₅, V, N)

    # Compute the weights for the higher-order cubature rule:
    WH₁ = (2 ^ N) * (12824 - 9120 * N + 400 * N ^ 2) / 19683
    WH₂ = (2 ^ N) * (980 / 6561)
    WH₃ = (2 ^ N) * (1820 - 400 * N) / 19683
    WH₄ = (2 ^ N) * (200 / 19683)
    WH₅ = 6859 / 19683

    # Compute the weights for the lower-order cubature rule:
    WL₁ = (2 ^ N) * (729 - 950 * N + 50 * N ^ 2) / 729
    WL₂ = (2 ^ N) * (245 / 486)
    WL₃ = (2 ^ N) * (265 - 100 * N) / 1458
    WL₄ = (2 ^ N) * (25 / 729)

    # Cache the results:
    GMCR         = GenzMalikCubatureRule{N}((P₂, P₃, P₄, P₅), (WH₁, WH₂, WH₃, WH₄, WH₅), (WL₁, WL₂, WL₃, WL₄))
    GMCRCache[N] = GMCR

    # Return the results:
    return GMCR
end

# Perform the integration using the Genz-Malik cubature rule:
function (GMCR::GenzMalikCubatureRule{N})(F, a::AbstractVector{T}, b::AbstractVector{T}) where {N, T<:Real}
    M = (b + a) / 2 
    Δ = (b - a) / 2 
    V = prod(Δ)

    # Evaluate the function at the center of the integration domain:
    F₁ = F(M)

    # Preallocate:
    F₂          = zero(F₁)
    F₃          = zero(F₁)
    F₄          = zero(F₁)
    F₅          = zero(F₁)
    EDimensions = similar(StaticArrays.SVector{N, typeof(LinearAlgebra.norm(F₁))})

    # Evaluate the function at the second and third set of integration points:
    for i = 1:N
        P₂  = Δ .* GMCR.P[1][i]
        F₂ᵢ = F(M + P₂) + F(M - P₂)

        P₃  = Δ .* GMCR.P[2][i]
        F₃ᵢ = F(M + P₃) + F(M - P₃)

        F₂  += F₂ᵢ
        F₃  += F₃ᵢ
        
        # Compute the error along each dimension:
        EDimensions[i] = LinearAlgebra.norm(F₃ᵢ + 12 * F₁ - 7*F₂ᵢ)
    end

    # Evaluate the function at the fourth of integration points:
    for P in GMCR.P[3]
        F₄ += F(M .+ Δ .* P)
    end

    # Evaluate the function at the fifth of integration points:
    for P in GMCR.P[4]
        F₅ += F(M .+ Δ .* P)
    end

    # Compute the integrals and error of integration:
    IH = V * (GMCR.WH[1] * F₁ + GMCR.WH[2] * F₂ + GMCR.WH[3] * F₃ + GMCR.WH[4] * F₄ + GMCR.WH[5] * F₅)
    IL = V * (GMCR.WL[1] * F₁ + GMCR.WL[2] * F₂ + GMCR.WL[3] * F₃ + GMCR.WL[4] * F₄)
    E  = LinearAlgebra.norm(IH - IL)

    # Find the dimension with the largest error:
    K = argmax(EDimensions)

    # Return the results:
    return IH, IL, E, K
end

function combinations(λ::T, ::Val{N}, K::Integer) where {N, T <: Real}
    Combinations = Combinatorics.combinations(1:N, K)

    P = Vector{StaticArrays.SVector{N, T}}(undef, length(Combinations))
    V = similar(StaticArrays.SVector{N, T})
    for (i, c) in enumerate(Combinations)
        V    .= 0
        V[c] .= λ
        P[i]  = V
    end

    # Return the results:
    return P
end

function signedcombinations(λ::T, ::Val{N}, K::Integer) where {N, T <: Real}
    Combinations = Combinatorics.combinations(1:N, K)
    Signs        = vec(collect.(Base.Iterators.product(Base.Iterators.repeated((-1, 1), K)...)))

    q = 2 ^ K
    P = Vector{StaticArrays.SVector{N, T}}(undef, length(Combinations) * q)
    V = similar(StaticArrays.SVector{N, T})
    for (i, c) in enumerate(Combinations)
        j     = (i - 1) * q + 1
        V    .= 0
        V[c] .= λ
        P[j]  = V

        # Use the gray code to flip one sign at a time:
        for s = 1:(q - 1)
            V        .= 0
            V[c]      = λ * Signs[s]
            P[j + s]  = V
        end
    end

    # Return the results:
    return P
end