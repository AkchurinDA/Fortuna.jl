# Sampling Random Variables

`Fortuna.jl` package allows to easily generate samples of both uncorrelated and correlated random variables using `samplerv()` function using different sampling techniques. Current version of the package implements [Inverse Transform Sampling (ITS)](https://en.wikipedia.org/wiki/Inverse_transform_sampling) and [Latin Hypercube Sampling (LHS)](https://en.wikipedia.org/wiki/Latin_hypercube_sampling) techniques.

```@setup 1
using Fortuna
using Random
Random.seed!(1)
```

## Sampling Random Variables

To generate samples of a random variable:

- Generate a random variable (`X`).

```@example 1
X = generaterv("Gamma", "M", [10, 1.5])

nothing # hide
```

- Sample the generated random variable using a sampling technique of your choice.

```@example 1
XSamples = samplerv(X, 5000, ITS())

nothing # hide
```

```@raw html
<img src="../Sample-RVariable.svg" class="center" style="width:600px; border-radius:5px;"/>
```

## Sampling Random Vectors with Uncorrelated Marginals

To generate samples of a random vector with *uncorrelated* marginals:

- Generate random variables (`X₁` and `X₂`).
- Define a random vector (`X`) with the generated random variables as marginals.

```@example 1
# Generate a random vector X with uncorrelated marginals X₁ and X₂:
X₁  = generaterv("Gamma", "M", [10, 1.5])
X₂  = generaterv("Gumbel", "M", [15, 2.5])
X   = [X₁, X₂]

nothing # hide
```

- Sample the defined random vector using a sampling technique of your choice.

```@example 1
# Generate 5000 samples of the random vector X using Inverse Transform Sampling technique:
XSamples = samplerv(X, 5000, ITS())

nothing # hide
```

```@raw html
<img src="../Sample-RVector-U.svg" class="center" style="width:450px; border-radius:5px;"/>
```

## Sampling Random Vectors with Correlated Marginals

To generate samples of a random vector with *correlated* marginals:

- Generate random variables (`X₁` and `X₂`).
- Define a random vector (`X`) with the generated random variables as marginals.

```@example 1
# Define a random vector X with correlated marginals X₁ and X₂:
X₁  = generaterv("Gamma", "M", [10, 1.5])
X₂  = generaterv("Gumbel", "M", [15, 2.5])
X   = [X₁, X₂]

nothing # hide
```

- Define a correlated matrix (`ρˣ`) for the defined random vector.
- Define a transformation object that hold all information about the define random vector.

```@example 1
# Define a correlation matrix:
ρˣ = [1 0.90; 0.90 1]

# Define a transformation object:
TransformationObject = NatafTransformation(X, ρˣ)

nothing # hide
```

Sample the defined random vector using a sampling technique of your choice.

```@example 1
# Generate 5000 samples of the random vector X in X-, Z-, and U-spaces using Inverse Transform Sampling technique:
XSamples, ZSamples, USamples = samplerv(TransformationObject, 5000, ITS())

nothing # hide
```

```@raw html
<img src="../Sample-RVector-C.svg" class="center" style="width:450px; border-radius:5px;"/>
```

## API

```@docs
ITS
LHS
samplerv(Samplers::Union{<:Distribution, Vector{<:Distribution}}, NumSamples::Integer, SamplingTechnique::AbstractSamplingTechnique)
samplerv(Object::NatafTransformation, NumSamples::Integer, SamplingTechnique::AbstractSamplingTechnique)
```