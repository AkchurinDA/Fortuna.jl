# Defining Reliability Problems

## Overview

In generally, 3 main "items" are always need to fully define a reliability problem and successfully analyze it to find the associated probability of failure ``P_{f}`` and reliability index ``\beta``:

- ``\vec{X}`` - Random vector with correlated non-normal marginals
- ``\rho^{X}`` - Correlation matrix
- ``g(\vec{X})`` - Limit state function

`Fortuna.jl` package uses these 3 "items" to fully define reliability problems using a custom `ReliabilityProblem()` type as shown in the example below.

```@setup 1
using Fortuna
```

```@example 1
# Generate a random vector X with correlated marginal random variables X₁ and X₂:
X₁  = randomvariable("Normal", "M", [10, 2])
X₂  = randomvariable("Normal", "M", [20, 5])
X   = [X₁, X₂]

# Define a correlation matrix for the random vector X:
ρˣ = [1 0.5; 0.5 1]

# Define a limit state function:
g(x) = x[1]^2 - 2 * x[2]

# Define a reliability problem using the provided information:
Problem = ReliabilityProblem(X, ρˣ, g)

nothing # hide
```

!!! note
    The definition of the limit state function ``g(\vec{X})`` in `Fortuna.jl` package only pertains to its form (e.g., whether it is linear, square, exponential, etc. in each variable). The information about the random variables involved in the reliability problem is carried in the random vector ``\vec{X}`` and its correlation matrix ``\rho^{X}``, that you use when defining a reliability problem using a custom `ReliabilityProblem()` type.

## Analyzing Reliability Problems

After defining the reliability problem, `Fortuna.jl` allows to easily solve it using a whole suite of First- and Second-Order Reliability Methods through a single `analyze()` function as shown in the example below.

```@example 1
# Solve the reliability problem using an imporved Hasofer-Lind-Rackwitz-Fiessler method:
Solution = solve(Problem, FORM(iHLRF()))
println("PoF = $(Solution.PoF)")
println("β = $(Solution.β)")
```

Descriptions of all First- and Second-Order Reliability Methods implemented in `Fortuna.jl` can be found on [First-Order Reliability Methods](@ref FORMPage) and [Second-Order Reliability Methods](@ref SORMPage) pages.

## API

```@docs
ReliabilityProblem
```