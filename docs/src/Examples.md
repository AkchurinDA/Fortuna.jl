# Examples

```@setup 1
using Fortuna
using Distributions
using Random
Random.seed!(123)
```

## Isoprobabilistic Transformations

### Nataf Transformation

```@example 1
# Define random vector:
X₁  = randomvariable("Gamma", "M", [10, 1.5])
X₂  = randomvariable("Gamma", "M", [15, 2.5])
X   = [X₁, X₂]

# Define correlation matrix:
ρˣ = [1 -0.75; -0.75 1]

# Perform Nataf Transformation:
NatafObject = NatafTransformation(X, ρˣ)

# Generate 10000 samples of random vector in X-, Z-, and U-spaces using Latin Hypercube Sampling technique:
XSamples, USamples, ZSamples = rand(NatafObject, 10000, LHS())

nothing # hide
```

```@raw html
<img src="../assets/Plots (Examples)/NatafTransformation-1.svg" class="center" style="max-height:350px; border-radius:2.5px;"/>
```

```@example 1
# Compute joint PDF of random vector:
xRange₁ = range(0, 20, 500)
xRange₂ = range(5, 25, 500)
fSamples = [pdf(NatafObject, [x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]

nothing # hide
```

```@raw html
<img src="../assets/Plots (Examples)/NatafTransformation-2.svg" class="center" style="max-height:350px; border-radius:2.5px;"/>
```

## Monte Carlo Methods

### Direct Monte Carlo Simulations

```@example 1
# Define random vector:
X₁  = randomvariable("Normal", "M", [0, 1])
X₂  = randomvariable("Normal", "M", [0, 1])
X   = [X₁, X₂]

# Define correlation matrix:
ρˣ  = [1 0; 0 1]

# Define limit state function:
β               = 3
g(x::Vector)    = β * sqrt(2) - x[1] - x[2]

# Define reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform reliability analysis using Monte Carlo simulations:
Solution = solve(Problem, MC())
println("MC:")
println("PoF: $(Solution.PoF)")
```

```@raw html
<img src="../assets/Plots (Examples)/MonteCarlo-1.svg" class="center" style="max-height:350px; border-radius:2.5px;"/>
```

### Importance Sampling Method

```@example 1
# Define proposal probability density function:
ProposalPDF = MvNormal([β / sqrt(2), β / sqrt(2)], [1 0; 0 1])

# Perform reliability analysis using Monte Carlo simulations:
Solution = solve(Problem, IS(q = ProposalPDF))
println("IS:")
println("PoF: $(Solution.PoF)")
```

```@raw html
<img src="../assets/Plots (Examples)/ImportanceSampling-1.svg" class="center" style="max-height:350px; border-radius:2.5px;"/>
```

```@raw html
<img src="../assets/Plots (Examples)/ImportanceSampling-2.svg" class="center" style="max-height:350px; border-radius:2.5px;"/>
```

## First-Order Reliability Methods

### Mean-Centered First-Order Second-Moment Method

```@example 1
# Define random vector:
X₁  = randomvariable("Normal", "M", [10, 2])
X₂  = randomvariable("Normal", "M", [20, 5])
X   = [X₁, X₂]

# Define correlation matrix:
ρˣ = [1 0.5; 0.5 1]

# Define two equivalent limit state functions:
g₁(x::Vector) = x[1] ^ 2 - 2 * x[2]
g₂(x::Vector) = 1 - 2 * x[2] / x[1] ^ 2

# Define reliability problems:
Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

# Perform reliability analysis using Mean-Centered First-Order Second-Moment (MCFOSM) method:
Solution₁ = solve(Problem₁, FORM(MCFOSM()))
Solution₂ = solve(Problem₂, FORM(MCFOSM()))
println("MCFOSM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")
```

### Hasofer-Lind-Rackwitz-Fiessler Method

```@example 1
# Perform reliability analysis using Hasofer-Lind-Rackwitz-Fiessler (HLRF) method:
Solution₁ = solve(Problem₁, FORM(HLRF()))
Solution₂ = solve(Problem₂, FORM(HLRF()))
println("FORM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")
```

### Improved Hasofer-Lind-Rackwitz-Fiessler Method

```@example 1
# Perform reliability analysis using improved Hasofer-Lind-Rackwitz-Fiessler (iHLRF) method:
Solution₁ = solve(Problem₁, FORM(iHLRF()))
Solution₂ = solve(Problem₂, FORM(iHLRF()))
println("FORM:")
println("β from g₁: $(Solution₁.β)")
println("β from g₂: $(Solution₂.β)")
```

## Second-Order Reliability Methods

### Curve-Fitting Method

```@example 1
# Define random vector:
M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]

# Define correlation matrix:
ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define limit state function:
a               = 0.190
s₁              = 0.030
s₂              = 0.015
g(x::Vector)    = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform reliability analysis using Curve-Fitting (CF) method:
Solution = solve(Problem, SORM(CF()))
println("SORM:")
println("β from FORM:       $(Solution.FORMSolution.β)")
println("β from SORM:       $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β from SORM:       $(Solution.β₂[2]) (Breitung)")
println("PoF from FORM:     $(Solution.FORMSolution.PoF)")
println("PoF from SORM:     $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF from SORM:     $(Solution.PoF₂[2]) (Breitung)")
```

### Point-Fitting Method

```@example 1
# Perform reliability analysis using point-fitting SORM:
Solution = solve(Problem, SORM(PF()))
println("SORM:")
println("β from FORM:       $(Solution.FORMSolution.β)")
println("β from SORM:       $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β from SORM:       $(Solution.β₂[2]) (Breitung)")
println("PoF from FORM:     $(Solution.FORMSolution.PoF)")
println("PoF from SORM:     $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF from SORM:     $(Solution.PoF₂[2]) (Breitung)")
```

## Subset Simulation Method

```@example 1
# Define random vector:
X₁  = randomvariable("Normal", "M", [0, 1])
X₂  = randomvariable("Normal", "M", [0, 1])
X   = [X₁, X₂]

# Define correlation matrix:
ρˣ  = [1 0; 0 1]

# Define limit state function:
β               = 3
g(x::Vector)    = β * sqrt(2) - x[1] - x[2]

# Define reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform reliability analysis using Monte Carlo simulations:
Solution = solve(Problem, SSM())
println("SSM:")
println("PoF: $(Solution.PoF)")
```

```@raw html
<img src="../assets/Plots (Examples)/SubsetSimulationMethod-1.svg" class="center" style="max-height:350px; border-radius:2.5px;"/>
```