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
Solution₁ = analyze(Problem₁, FORM(MCFOSM()))
Solution₂ = analyze(Problem₂, FORM(MCFOSM()))
println("MCFOSM:")
println("β from G₁: $(Solution₁.β)")
println("β from G₂: $(Solution₂.β)")

# Perform the reliability analysis using FORM:
Solution₁ = analyze(Problem₁, FORM(iHLRF()))
Solution₂ = analyze(Problem₂, FORM(iHLRF()))
println("FORM:")
println("β from G₁: $(Solution₁.β)")
println("β from G₂: $(Solution₂.β)")