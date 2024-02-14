# Examples

```@setup 1
using Fortuna
using Random
Random.seed!(1)
```

## Example #1: Perform Nataf Transformation

```@example 1
# Generate a random vector X with correlated marginal random variables X₁ and X₂:
X₁  = generaterv("Gamma", "M", [10, 1.5])
X₂  = generaterv("Gumbel", "M", [15, 2.5])
X   = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0.90; 0.90 1]

# Perform the Nataf Transformation by defining a "NatafTransformation" object:
NatafObject = NatafTransformation(X, ρˣ)

# Generate 1000 samples of the random vector X in X-, Z-, and U-spaces:
XSamples, USamples, ZSamples = samplerv(NatafObject, 1000, ITS())

nothing # hide
```

![Nataf Transformation](./assets/NatafTransformation.svg)

## Example #3: Reliability Analysis (FORM - MCFOSM)

```@example 1
# Generate a random vector X with correlated marginal random variables:
X₁  = generaterv("Normal", "Moments", [10, 2])
X₂  = generaterv("Normal", "Moments", [20, 5])
X   = [X₁, X₂]
ρˣ  = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
g₁(x) = x[1] ^ 2 - 2 * x[2]
g₂(x) = 1 - 2 * x[2] / x[1] ^ 2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

# Perform the reliability analysis:
Solution₁ = analyze(Problem₁, FORM(MCFOSM()))
Solution₂ = analyze(Problem₂, FORM(MCFOSM()))
println("MCFOSM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")
```

## Example #4: Reliability Analysis (FORM - iHLRF)

```@example 1
# Generate a random vector X with correlated marginal random variables:
X₁  = generaterv("Normal", "Moments", [10, 2])
X₂  = generaterv("Normal", "Moments", [20, 5])
X   = [X₁, X₂]
ρˣ  = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
g₁(x) = x[1] ^ 2 - 2 * x[2]
g₂(x) = 1 - 2 * x[2] / x[1] ^ 2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

# Perform the reliability analysis:
Solution₁ = analyze(Problem₁, FORM(iHLRF()))
Solution₂ = analyze(Problem₂, FORM(iHLRF()))
println("FORM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")

```

## Example #5: Reliability Analysis (MCS)

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁  = generaterv("Normal", "M", [250, 250 * 0.3])
M₂  = generaterv("Normal", "M", [125, 125 * 0.3])
P   = generaterv("Gumbel", "M", [2500, 2500 * 0.2])
Y   = generaterv("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a   = 0.190
s₁  = 0.030
s₂  = 0.015
g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using curve-fitting SORM:
Solution = analyze(Problem, MCS())
println("MCS:")
println("PoF: $(Solution.PoF)")
```

## Example #5: Reliability Analysis (SORM - CF)

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁  = generaterv("Normal", "M", [250, 250 * 0.3])
M₂  = generaterv("Normal", "M", [125, 125 * 0.3])
P   = generaterv("Gumbel", "M", [2500, 2500 * 0.2])
Y   = generaterv("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a   = 0.190
s₁  = 0.030
s₂  = 0.015
g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

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

## Example #6: Reliability Analysis (SORM - PF)

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁  = generaterv("Normal", "M", [250, 250 * 0.3])
M₂  = generaterv("Normal", "M", [125, 125 * 0.3])
P   = generaterv("Gumbel", "M", [2500, 2500 * 0.2])
Y   = generaterv("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a   = 0.190
s₁  = 0.030
s₂  = 0.015
g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using curve-fitting SORM:
Solution = analyze(Problem, SORM(PF()))
println("SORM:")
println("β from FORM: $(Solution.β₁)")
println("β from SORM: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β from SORM: $(Solution.β₂[2]) (Breitung)")
println("PoF from FORM: $(Solution.PoF₁)")
println("PoF from SORM: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF from SORM: $(Solution.PoF₂[2]) (Breitung)")
```

## Example #7: Reliability Analysis (SSM)
```@example 1
# Define a random vector of uncorrelated marginal distributions:
X₁  = generaterv("Exponential", "P", 1)
X₂  = generaterv("Exponential", "P", 1)
X   = [X₁, X₂]
ρˣ  = [1 0; 0 1]

# Define a limit state function:
g(x) = 10 - x[1] - x[2]

# Define reliability problems:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using SSM:
Solution = analyze(Problem, SSM())
println("SSM")
println("PoF from SSM: $(Solution.PoF)")
```