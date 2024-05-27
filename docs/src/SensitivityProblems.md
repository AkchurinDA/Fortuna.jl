# Sensitivity Problems

## Types of Sensitivity Problems in `Fortuna.jl`

`Fortuna.jl` districts two types of sensitivity problems:

| Type | Description |
| :--- | :--- |
| I | Used to find sensitivities w.r.t. to the parameters ``\vec{\Theta}_{g}`` of the limit state function ``g(\vec{X}, \vec{\Theta}_{g})`` |
| II | Used to find sensitivities w.r.t. to the parameters and/or moments ``\vec{\Theta}_{f}`` of the random vector ``\vec{X}(\vec{\Theta}_{f})`` |

## Defining and Solving Sensitivity Problems of Type I

In general, 4 main "items" are always needed to fully define a sensitivity problem of type I and successfully solve it to find the associated sensitivity vectors of probability of failure ``\vec{\nabla}_{\vec{\Theta}_{g}} P_{f}`` and reliability index ``\vec{\nabla}_{\vec{\Theta}_{g}} \beta`` w.r.t. to the parameters of the limit state function ``\vec{\Theta}_{g}``:

| Item | Description |
| :--- | :--- |
| ``\vec{X}`` | Random vector |
| ``\rho^{X}`` | Correlation matrix |
| ``g(\vec{X}, \vec{\Theta}_{g})`` | Limit state function parametrized in terms of its parameters |
| ``\vec{\Theta}_{g}`` | Parameters of the limit state function |

`Fortuna.jl` package uses these 4 "items" to fully define sensitivity problems of type I using `SensitivityProblem()` type as shown in the example below.

```@setup 1
using Fortuna
```

```@example 1
# Define the random vector:
M₁ = randomvariable("Normal",  "M", [250,   250   * 0.3])
M₂ = randomvariable("Normal",  "M", [125,   125   * 0.3])
P  = randomvariable("Gumbel",  "M", [2500,  2500  * 0.2])
Y  = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
X  = [M₁, M₂, P, Y]

# Define the correlation matrix:
ρˣ = [
    1.0 0.5 0.3 0.0
    0.5 1.0 0.3 0.0
    0.3 0.3 1.0 0.0
    0.0 0.0 0.0 1.0]

# Define the limit state function:
g(x::Vector, θ::Vector) = 1 - x[1] / (θ[1] * x[4]) - x[2] / (θ[2] * x[4]) - (x[3] / (θ[3] * x[4])) ^ 2

# Define parameters of the limit state function:
s₁ = 0.030
s₂ = 0.015
a  = 0.190
Θ  = [s₁, s₂, a]

# Define a sensitivity problem:
Problem = SensitivityProblemTypeI(X, ρˣ, g, Θ)

nothing # hide
```

After defining a sensitivity problem of type I, `Fortuna.jl` allows to easily perform sensitivity analysis using a single `solve()` function as shown in the example below.

```@example 1
# Perform the sensitivity analysis:
Solution = solve(Problem)
println("∇β   = $(Solution.∇β)  ")
println("∇PoF = $(Solution.∇PoF)")
```

## Defining and Solving Sensitivity Problems of Type II

Similar to sensitivity problem of type I, 4 main "items" are needed to fully define a sensitivity problem of type II and successfully solve it to find the associated sensitivity vectors of probability of failure ``\vec{\nabla}_{\vec{\Theta}_{f}} P_{f}`` and reliability index ``\vec{\nabla}_{\vec{\Theta}_{f}} \beta`` w.r.t. to the parameters and/or moments of the random vector ``\vec{\Theta}_{f}``:

| Item | Description |
| :--- | :--- |
| ``\vec{X}(\vec{\Theta}_{f})`` | Random vector with correlated non-normal marginals parameterized in terms of its parameters and/or moments |
| ``\rho^{X}`` | Correlation matrix |
| ``g(\vec{X})`` | Limit state function |
| ``\vec{\Theta}_{f}`` | Parameters and/or moments of the random vector |

`Fortuna.jl` package uses these 4 "items" to fully define sensitivity problems of type I using `SensitivityProblem()` type as shown in the example below.

```@setup 1
using Fortuna
```

```@example 1
# Define the random vector as a function of its parameters and moments:
function XFunction(Θ::Vector)
    M₁ = randomvariable("Normal",  "M", [Θ[1], Θ[2]])
    M₂ = randomvariable("Normal",  "M", [Θ[3], Θ[4]])
    P  = randomvariable("Gumbel",  "M", [Θ[5], Θ[6]])
    Y  = randomvariable("Weibull", "M", [Θ[7], Θ[8]])

    return [M₁, M₂, P, Y]
end

# Define the correlation matrix:
ρˣ = [
    1.0 0.5 0.3 0.0
    0.5 1.0 0.3 0.0
    0.3 0.3 1.0 0.0
    0.0 0.0 0.0 1.0]

# Define the parameters and moments of the random vector:
Θ = [
      250,   250 * 0.30,
      125,   125 * 0.30,
     2500,  2500 * 0.20,
    40000, 40000 * 0.10]

# Define the limit state function:
a            = 0.190
s₁           = 0.030
s₂           = 0.015
g(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

# Define a sensitivity problem:
Problem  = SensitivityProblemTypeII(XFunction, ρˣ, g, Θ)

nothing # hide
```

Similar to sensitivity problems of type I, sensitivity problems of type II are solved using the same `solve()` function as shown in the example below.

```@example 1
# Perform the sensitivity analysis:
Solution = solve(Problem)
println("∇β   = $(Solution.∇β)  ")
println("∇PoF = $(Solution.∇PoF)")
```

## API

```@docs
solve(Problem::SensitivityProblemTypeI)
solve(Problem::SensitivityProblemTypeII)
SensitivityProblemTypeI
SensitivityProblemTypeII
SensitivityProblemCache
```