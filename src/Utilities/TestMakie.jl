using Fortuna
using CairoMakie
CairoMakie.activate!(type = :svg)

# Define a random vector of correlated marginal distributions:
X₁ = randomvariable("Normal", "M", [200, 20])
X₂ = randomvariable("Normal", "M", [150, 10])
X  = [X₁, X₂]
ρˣ = [1 0; 0 1]

# Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
g(x::Vector) = x[1] - x[2]

# Define reliability problems:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using HLRF:
Solution = solve(Problem, FORM())

# Plot:
reliabilityplot(Problem, Solution)

typeof(Problem)
typeof(Solution)