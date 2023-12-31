# Generating Random Variables

`Fortuna.jl` package builds its capacity to generate random variables using `generaterv()` function by utilizing the widely-adopted [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) package, enabling seamless integration with other probabilistic programming Julia packages such as [`Turing.jl`](https://github.com/TuringLang/Turing.jl) and [`RxInfer.jl`](https://github.com/biaslab/RxInfer.jl). However, unlike [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) package, `Fortuna.jl` allows you to generate random variables not only using their *parameters*, but also using their *moments*, which often useful in the field of Structural and System Reliability Analysis.

```@docs
generaterv(DistributionName::String, DefineBy::String, Values::Union{Real,Vector{<:Real}})
```

```@setup GenerateRV
using Fortuna
```

## Generating Random Variables Using Moments

```@example GenerateRV
# Generate a lognormally distributed random variable R 
# with mean (μ) of 15 and standard deviation (σ) of 10:
R = generaterv("Lognormal", "M", [15, 10])
println("μ = $(mean(R))")
println("σ = $(std(R))")
```

## Generating Random Variables Using Parameters

```@example GenerateRV
# Generate a gamma-distributed random variable Q 
# with shape parameter (α) of 16 and scale parameter (θ) of 0.625:
Q = generaterv("Gamma", "P", [16, 0.625])
println("α = $(params(Q)[1])")
println("θ = $(params(Q)[2])")
```

## Supported Random Variables

`Fortuna.jl` package currently supports the following distributions:
- Exponential
- Gamma
- Gumbel
- Lognormal
- Normal
- Poisson
- Uniform
- Weibull

!!! tip
    If you want to define a random variable that is not supported by Fortuna.jl package, please raise an issue on the Github Issues page.