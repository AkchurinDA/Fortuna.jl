# Examples

```@setup 1
using Fortuna
using Random
Random.seed!(1)
```

## Example #1: Perform Nataf Transformation

```@example 1
# Generate a random vector X with correlated marginal random variables X₁ and X₂:
X₁ = generaterv("Gamma", "M", [10, 1.5])
X₂ = generaterv("Gumbel", "M", [15, 2.5])
X = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0.90; 0.90 1]

# Perform the Nataf Transformation by defining a "NatafTransformation" object:
NatafObject = NatafTransformation(X, ρˣ)

# Generate 1000 samples of the random vector X in X-, Z-, and U-spaces:
XSamples, USamples, ZSamples = samplerv(NatafObject, 1000)

nothing # hide
```

![Nataf Transformation](./assets/NatafTransformation.svg)

## Example #2: Reliability Analysis (FORM - MCFOSM)

```@example 1
# Generate a random vector X with correlated marginal random variables:
X₁ = generaterv("Normal", "Moments", [10, 2])
X₂ = generaterv("Normal", "Moments", [20, 5])
X = [X₁, X₂]
ρˣ = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
G₁(x::Vector) = x[1]^2 - 2 * x[2]
G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

# Perform the reliability analysis:
Solution₁ = analyze(Problem₁, FORM(MCFOSM()))
Solution₂ = analyze(Problem₂, FORM(MCFOSM()))
println("MCFOSM:")
println("β from G₁: $(Solution₁.β)")
println("β from G₂: $(Solution₂.β)")
```

## Example #3: Reliability Analysis (FORM - iHLRF)

```@example 1
# Generate a random vector X with correlated marginal random variables:
X₁ = generaterv("Normal", "Moments", [10, 2])
X₂ = generaterv("Normal", "Moments", [20, 5])
X = [X₁, X₂]
ρˣ = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
G₁(x::Vector) = x[1]^2 - 2 * x[2]
G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

# Perform the reliability analysis:
Solution₁ = analyze(Problem₁, FORM(iHLRF()))
Solution₂ = analyze(Problem₂, FORM(iHLRF()))
println("FORM:")
println("β from G₁: $(Solution₁.β)")
println("β from G₂: $(Solution₂.β)")
```

## Example #4: Reliability Analysis (SORM - CF)

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁ = generaterv("Normal", "Moments", [250, 250 * 0.3])
M₂ = generaterv("Normal", "Moments", [125, 125 * 0.3])
P = generaterv("Gumbel", "Moments", [2500, 2500 * 0.2])
Y = generaterv("Weibull", "Moments", [40000, 40000 * 0.1])
X = [M₁, M₂, P, Y]
ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a = 0.190
s₁ = 0.030
s₂ = 0.015
G(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, G)

# Perform the reliability analysis using curve-fitting SORM:
Solution = analyze(Problem, SORM(CF()))
println("SORM:")
println("β from FORM: $(Solution.β₁)")
println("β from SORM: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β from SORM: $(Solution.β₂[2]) (Breitung)")
println("PoF from FORM: $(Solution.PoF₁)")
println("PoF from SORM: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF from SORM: $(Solution.PoF₂[2]) (Breitung)")
```