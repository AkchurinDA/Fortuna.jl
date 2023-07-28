using Fortuna

# Define a random vector of correlated marginal distributions:
X₁ = generaterv("Normal", "Moments", [10, 2])
X₂ = generaterv("Normal", "Moments", [20, 5])
X = [X₁, X₂]
ρˣ = [1 0.5; 0.5 1]

# Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
G₁(x::Vector) = x[1]^2 - 2 * x[2]
G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

# Perform the reliability analysis using MCFOSM:
β₁ = analyze(Problem₁, MCFOSM())
β₂ = analyze(Problem₂, MCFOSM())
println("MCFOSM:")
println("β from G₁: $β₁")
println("β from G₂: $β₂")

# Perform the reliability analysis using FORM:
β₁, x₁, u₁ = analyze(Problem₁, FORM(iHLRF()))
β₂, x₂, u₂ = analyze(Problem₂, FORM(iHLRF()))
println("FORM:")
println("β from G₁: $β₁")
println("β from G₂: $β₂")