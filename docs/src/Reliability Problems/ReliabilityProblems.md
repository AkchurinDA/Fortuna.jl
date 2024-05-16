# Reliability Problems

## Defining Reliability Problems

In general, 3 main "items" are always need to fully define a reliability problem and successfully solve it to find the associated probability of failure ``P_{f}`` and reliability index ``\beta``:

| Item | Description |
| :--- | :--- |
| ``\vec{X}`` | Random vector with correlated non-normal marginals |
| ``\rho^{X}`` | Correlation matrix |
| ``g(\vec{X})`` | Limit state function |

`Fortuna.jl` package uses these 3 "items" to fully define reliability problems using a custom `ReliabilityProblem()` type as shown in the example below.

```@setup 1
using Fortuna
```

```@example 1
# Define random vector:
X₁  = randomvariable("Normal", "M", [10, 2])
X₂  = randomvariable("Normal", "M", [20, 5])
X   = [X₁, X₂]

# Define correlation matrix:
ρˣ = [1 0.5; 0.5 1]

# Define limit state function:
g(x::Vector) = x[1] ^ 2 - 2 * x[2]

# Define reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

nothing # hide
```

!!! note
    The definition of the limit state function ``g(\vec{X})`` in `Fortuna.jl` package only pertains to its form (e.g., whether it is linear, square, exponential, etc. in each variable). The information about the random variables involved in the reliability problem is carried in the random vector ``\vec{X}`` and its correlation matrix ``\rho^{X}``, that you use when defining a reliability problem using a custom `ReliabilityProblem()` type.

## Solving Reliability Problems

After defining the reliability problem, `Fortuna.jl` allows to easily solve it using a whole suite of First- and Second-Order Reliability Methods through a single `solve()` function as shown in the example below.

```@example 1
# Perform reliability analysis using improved Hasofer-Lind-Rackwitz-Fiessler (iHLRF) method:
Solution = solve(Problem, FORM(iHLRF()))
println("β   = $(Solution.β)  ")
println("PoF = $(Solution.PoF)")
```

Descriptions of all First- and Second-Order Reliability Methods implemented in `Fortuna.jl` can be found on [First-Order Reliability Methods](@ref FORMPage) and [Second-Order Reliability Methods](@ref SORMPage) pages.

## API

```@docs
ReliabilityProblem
```