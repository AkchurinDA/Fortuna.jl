# Inverse Reliability Problems

## Defining Inverse Reliability Problems

In general, 4 main "items" are always need to fully define an inverse reliability problem and successfully solve it to find the parameter of interest ``\theta``:

| Item | Description |
| :--- | :--- |
| ``\vec{X}`` | Random vector with correlated non-normal marginals |
| ``\rho^{X}`` | Correlation matrix |
| ``g(\vec{X}, \theta)`` | Limit state function |
| ``\beta`` | Target reliability index |

`Fortuna.jl` package uses these 4 "items" to fully define inverse reliability problems of type I using `SensitivityProblem()` type as shown in the example below.

## Solving Inverse Reliability Problems

```@setup InverseReliabilityProblem
using Fortuna
```

```@example InverseReliabilityProblem
# Define the random vector:
X₁ = randomvariable("Normal", "M", [0, 1])
X₂ = randomvariable("Normal", "M", [0, 1])
X₃ = randomvariable("Normal", "M", [0, 1])
X₄ = randomvariable("Normal", "M", [0, 1])
X  = [X₁, X₂, X₃, X₄]

# Define the correlation matrix:
ρˣ = Matrix{Float64}(1.0 * I, 4, 4)

# Define the limit state function:
g(x::Vector, θ::Real) = exp(-θ * (x[1] + 2 * x[2] + 3 * x[3])) - x[4] + 1.5

# Define the target reliability index:
β = 2

# Define an inverse reliability problem:
Problem = InverseReliabilityProblem(X, ρˣ, g, β)

nothing # hide
```

## Solving Inverse Reliability Problems

After defining an inverse reliability problem, `Fortuna.jl` allows to easily solve it using a single `solve()` function as shown in the example below.

```@example InverseReliabilityProblem
# Perform the inverse reliability analysis:
Solution = solve(Problem, 0.1, x₀ = [0.2, 0.2, 0.2, 0.2])
println("x = $(Solution.x[:, end])  ")
println("θ = $(Solution.θ[end])")
```

## API

```@docs
solve(Problem::InverseReliabilityProblem, θ₀::Real)
InverseReliabilityProblem
InverseReliabilityProblemCache
```