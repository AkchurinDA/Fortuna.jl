# Examples

```@setup 1
using Fortuna
using Random
Random.seed!(1)
```

## Isoprobabilistic Transformation

### Nataf Transformation

```@example 1
# Generate a random vector X with correlated marginal random variables X₁ and X₂:
X₁  = randomvariable("Gamma", "M", [10, 1.5])
X₂  = randomvariable("Gumbel", "M", [15, 2.5])
X   = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0.90; 0.90 1]

# Perform the Nataf Transformation by defining a "NatafTransformation" object:
NatafObject = NatafTransformation(X, ρˣ)

# Generate 1000 samples of the random vector X in X-, Z-, and U-spaces:
XSamples, USamples, ZSamples = rand(NatafObject, 1000, ITS())

nothing # hide
```

```@raw html
<img src="../assets/NatafTransformation.svg" class="center" style="border-radius:5px;"/>
```

## Monte Carlo Simulations

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a   = 0.190
s₁  = 0.030
s₂  = 0.015
g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using Monte Carlo simulations:
Solution = solve(Problem, MC())
println("MC:")
println("PoF: $(Solution.PoF)")
```

## First-Order Reliability Methods

### Mean-Centered First-Order Second-Moment Method

```@example 1
# Generate a random vector X with correlated marginal random variables:
X₁  = randomvariable("Normal", "M", [10, 2])
X₂  = randomvariable("Normal", "M", [20, 5])
X   = [X₁, X₂]
ρˣ  = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
g₁(x) = x[1] ^ 2 - 2 * x[2]
g₂(x) = 1 - 2 * x[2] / x[1] ^ 2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

# Perform the reliability analysis:
Solution₁ = solve(Problem₁, FORM(MCFOSM()))
Solution₂ = solve(Problem₂, FORM(MCFOSM()))
println("MCFOSM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")
```

### Improved Hasofer-Lind-Rackwitz-Fiessler Method

```@example 1
# Generate a random vector X with correlated marginal random variables:
X₁  = randomvariable("Normal", "M", [10, 2])
X₂  = randomvariable("Normal", "M", [20, 5])
X   = [X₁, X₂]
ρˣ  = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
g₁(x) = x[1] ^ 2 - 2 * x[2]
g₂(x) = 1 - 2 * x[2] / x[1] ^ 2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

# Perform the reliability analysis:
Solution₁ = solve(Problem₁, FORM(iHLRF()))
Solution₂ = solve(Problem₂, FORM(iHLRF()))
println("FORM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")
```

## Second-Order Reliability Methods

### Curve-Fitting Method

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a   = 0.190
s₁  = 0.030
s₂  = 0.015
g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using Curve-Fitting (CF) method:
Solution = solve(Problem, SORM(CF()))
println("SORM:")
println("β from FORM: $(Solution.FORMSolution.β)")
println("β from SORM: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β from SORM: $(Solution.β₂[2]) (Breitung)")
println("PoF from FORM: $(Solution.FORMSolution.PoF)")
println("PoF from SORM: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF from SORM: $(Solution.PoF₂[2]) (Breitung)")
```

### Point-Fitting Method

```@example 1
# Generate a random vector X with correlated marginal random variables:
M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define a limit state function:
a   = 0.190
s₁  = 0.030
s₂  = 0.015
g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using Point-Fitting (PF) method:
Solution = solve(Problem, SORM(PF()))
println("SORM:")
println("β from FORM: $(Solution.FORMSolution.β)")
println("β from SORM: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β from SORM: $(Solution.β₂[2]) (Breitung)")
println("PoF from FORM: $(Solution.FORMSolution.PoF)")
println("PoF from SORM: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF from SORM: $(Solution.PoF₂[2]) (Breitung)")
```