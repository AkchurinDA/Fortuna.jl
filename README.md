# STRELA
[![Build Status](https://github.com/AkchurinDA/STRELA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AkchurinDA/STRELA.jl/actions/workflows/CI.yml)

## Description
`STRELA` (**S**tructural **Rel**iability **A**nalysis) is a general purpose Julia package for structural and system reliability analysis.

<p align="center">
  <img src="assets/Logo.svg" alt = "Logo" width="25%">
</p>

## Installation
To install `STRELA` package, type `]` in Julia REPL to enter package manager mode and execute the following command:

```
pkg> add STRELA
```

## Quick Start
To start using `STRELA` package, type the following command in Julia REPL or in the beginning of a file:

```julia
using STRELA
```

### Generating Random Variables
`STRELA` package builds its capacity to generate random variables using `generaterv()` function by utilizing the widely-adopted [Distributions](https://github.com/JuliaStats/Distributions.jl) package, enabling seamless integration with other Julia packages such as [Turing](https://github.com/TuringLang/Turing.jl). However, unlike [Distributions](https://github.com/JuliaStats/Distributions.jl) package, STRELA allows you to generate random variables not only using their **parameters**, but also using their **moments**, which often useful.

```julia
# Generate a lognormally distributed random variable R with mean (μ) of 15 and standard deviation (σ) of 2.5:
R = generaterv("Lognormal", "Moments", [15, 2.5])

# Generate a gamma-distributed random variable Q with shape parameter (α) of 16 and scale parameter (θ) of 0.625:
Q = generaterv("Gamma", "Parameters", [16, 0.625])
```

### Performing Nataf Transformation
`STRELA` package allows to easily perform the Nataf transformation of correlated random variables into the space of uncorrelated standard normal variables.

```julia
# Define a random vector:
X₁ = generaterv("Gamma", "Moments", [10, 1.5])
X₂ = generaterv("Gumbel", "Moments", [15, 2.5])
X = [X₁, X₂]

# Define correlation coefficients between marginal random variables of the random vector:
ρˣ = [1 0.75; 0.75 1]

# Perform Nataf transformation of the random vector by defining a "NatafTransformation" object:
NatafObject = NatafTransformation(X, ρˣ)
```

The results of the performed Nataf transformation can be accesses from the fields of the defined `NatafTransformation` object.

```julia
# Extract the distorted correlation matrix of correlated normal random variables:
display(NatafObject.ρᶻ)
# 2×2 Matrix{Float64}:
#  1.0       0.765315
#  0.765315  1.0

# Extract the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix:
display(NatafObject.L)
# 2×2 LinearAlgebra.LowerTriangular{Float64, Matrix{Float64}}:
#  1.0        ⋅ 
#  0.765315  0.643656

# Extract the inverse of the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix
display(NatafObject.L⁻¹)
# 2×2 LinearAlgebra.LowerTriangular{Float64, Matrix{Float64}}:
#   1.0       ⋅ 
#  -1.18901  1.55362
```

### Sampling Random Variables
`STRELA` package also allows to easily generate samples of uncorrelated and correlated random variables using `samplerv()` function using different sampling techniques. Current version of the package implements [Inverse Transform Sampling (ITS)](https://en.wikipedia.org/wiki/Inverse_transform_sampling) and [Latin Hypercube Sampling (LHS)](https://en.wikipedia.org/wiki/Latin_hypercube_sampling) techniques.

#### Uncorrelated Random Variables
The function `generaterv()` allows to generate samples of a single distribution, as well as to generate samples of random vectors.
```julia
# Define a random vector:
X₁ = generaterv("Gamma", "Moments", [10, 1.5])
X₂ = generaterv("Gumbel", "Moments", [15, 2.5])
X = [X₁, X₂]

# Generate 3 samples of the random variable X₁ using Inverse Transform Sampling:
X₁SamplesITS = samplerv(X₁, 3, "ITS")
# 3-element Vector{Float64}:
#  8.438113227625095
#  9.103174415760643
#  11.11171748034975

# Generate 3 samples of the random variable X₁ using Latin Hypercube Sampling:
X₁SamplesLHS = samplerv(X₁, 3, "LHS")
# 3-element Vector{Float64}:
#  10.70297332783710
#  9.380731608864231
#  8.997067439248992

# Generate 3 samples of the random vector using Inverse Transform Sampling:
XSamplesITS = samplerv(X, 3, "ITS")
# 3×2 Matrix{Float64}:
#  11.1815  11.6162
#  11.0042  12.6362
#  10.4576  13.8437

# Generate 3 samples of the random vector using Latin Hypercube Sampling:
XSamplesLHS = samplerv(X, 3, "LHS")
# 3×2 Matrix{Float64}:
#  9.14107  12.6207
#  9.64263  13.5258
#  8.63828  19.3058
```

#### Correlated Random Variables
Generating the correlated random variables can be done by:
1. Performing the Nataf transformation of the random correlated variables.

    ```julia
    # Define a random vector:
    X₁ = generaterv("Gamma", "Moments", [10, 1.5])
    X₂ = generaterv("Gumbel", "Moments", [15, 2.5])
    X = [X₁, X₂]

    # Define correlation coefficients between marginal distributions of the random vector:
    ρˣ = [1 0.75; 0.75 1]

    # Perform Nataf transformation by defining a "NatafTransformation" object:
    NatafObject = NatafTransformation(X, ρˣ)
    ```

2. Generating samples of the random vector with correlated marginal random variables by passing the defined `NatafTransformation` object directly into the sampling function `samplerv()`.

    ```julia
    # Generate 3 samples of the random vector using Inverse Transform Sampling:
    XSamplesITS, _, _ = samplerv(NatafObject, 3, "ITS")
    #  10.1631  14.1579
    #  9.42206  12.6114
    #  12.3663  14.9653

    # Generate 3 samples of the random vector using Latin Hypercube Sampling:
    XSamplesLHS, _, _ = samplerv(NatafObject, 3, "LHS")
    # 3×2 Matrix{Float64}:
    #  11.3783  17.6984
    #  9.02893  12.6861
    #  10.0632  14.2938
    ```

### Reliability Analysis 
Ultimately, `STRELA` package is developed to perform structural and system reliability analysis. The current version of the package implements the Mean-Centered First-Order Second-Moment (MCFOSM) reliability method, as well as the First-Order Reliability Method (FORM).

#### Mean-Centered First-Order Second-Moment (MCFOSM) Reliability Method
The MCFOSM method is the simplest and least expensive type of reliability method. It utilizes the first-order Taylor expansion of the limit state function at the mean values and the first two moments of the random variables involved in the reliability problem to evaluate the reliability index. However, despite the fact that it is simple and does not require the complete knowledge of the random variables involved in the reliability problem, the MCFOSM method faces an issue known as the invariance problem. This problem arises because the resulting reliability index is dependent on the formulation of the limit state function. In other words, two equivalent limit state functions with the same failure boundaries produce two different reliability indices.

```julia
# Define a random vector of correlated marginal distributions:
X₁ = generaterv("Normal", "Moments", [10, 2])
X₂ = generaterv("Normal", "Moments", [20, 5])
X = [X₁, X₂]
ρˣ = [1 0.5; 0.5 1]

# Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
G₁(x::Vector) = x[1]^2 - 2 * x[2]
G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

# Perform the reliability analysis using MCFOSM:
β₁ = MCFOSM(Problem₁)
β₂ = MCFOSM(Problem₂)
println("MCFOSM:")
println("β from G₁: $β₁")
println("β from G₂: $β₂")
# MCFOSM:
# β from G₁: 1.664100588675687
# β from G₂: 4.285714285714286
```

#### First-Order Reliability Method (FORM)
The FORM overcomes the invariance problem faced by the MCFOSM method by using the first-order Taylor expansion of the limit state function at a point known as the "design point" on the failure boundary. Since the design point is not known a priori, the FORM is inherently an iterative method. Current implementation of FORM uses improved Hasofer-Lind Rackwitz-Fiessler (iHL-RF) algorithm.

```julia
β₁, _, _ = FORM(Problem₁)
β₂, _, _ = FORM(Problem₂)
println("FORM:")
println("β from G₁: $β₁")
println("β from G₂: $β₂")
# FORM:
# β from G₁: 2.10833940741697
# β from G₂: 2.10833972384163
```

# License
`STRELA` package is distributed under the [MIT license](https://en.wikipedia.org/wiki/MIT_License). More information can be found in the `LICENSE` file.

# Help and Support
For assistance with the package, please raise an issue on the Github Issues page. Please use the appropriate labels to indicate the specific functionality you are inquiring about.