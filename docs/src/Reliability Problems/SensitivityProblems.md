# Sensitivity Problems

## Defining Sensitivity Problems

In general, 4 main "items" are always need to fully define a sensitivity problem and successfully solve it to find the associated sensitivity vectors of probability of failure ``\vec{\nabla}_{\vec{\Theta}} P_{f}`` and reliability index ``\vec{\nabla}_{\vec{\Theta}} \beta`` with respect to limit state function's parameters ``\vec{\Theta}``:

| Item | Description |
| :--- | :--- |
| ``\vec{X}`` | Random vector with correlated non-normal marginals |
| ``\rho^{X}`` | Correlation matrix |
| ``g(\vec{X}, \vec{\Theta})`` | Limit state function |
| ``\vec{\theta}`` | Parameters of limit state function |

`Fortuna.jl` package uses these 4 "items" to fully define sensitivity problems using a custom `SensitivityProblem()` type as shown in the example below.

```@setup 1
using Fortuna
```

```@example 1
# Define random vector:
M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
X   = [M₁, M₂, P, Y]

# Define correlation matrix:
ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

# Define limit state function:
g(x::Vector, θ::Vector) = 1 - x[1] / (θ[1] * x[4]) - x[2] / (θ[2] * x[4]) - (x[3] / (θ[3] * x[4])) ^ 2

# Define parameters of limit state function:
s₁  = 0.030
s₂  = 0.015
a   = 0.190
θ   = [s₁, s₂, a]

# Define sensitivity problem:
Problem = SensitivityProblem(X, ρˣ, g, θ)

nothing # hide
```

## Solving Sensitivity Problems

After defining the sensitivity problem, `Fortuna.jl` allows to easily solve it using a single `solve()` function as shown in the example below.

```@example 1
# Perform reliability analysis using improved Hasofer-Lind-Rackwitz-Fiessler (iHLRF) method:
Solution = solve(Problem)
println("∇PoF   = $(Solution.∇PoF)")
println("∇β     = $(Solution.∇β)")
```

## API

```@docs
solve(Problem::SensitivityProblem)
SensitivityProblem
SensitivityProblemCache
```