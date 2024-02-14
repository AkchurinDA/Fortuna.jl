# Generating Random Variables

`Fortuna.jl` package builds its capacity to generate random variables using `generaterv()` function by utilizing the widely-adopted [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) package, enabling seamless integration with other probabilistic programming Julia packages such as [`Turing.jl`](https://github.com/TuringLang/Turing.jl) and [`RxInfer.jl`](https://github.com/biaslab/RxInfer.jl). However, unlike [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) package, `Fortuna.jl` allows you to generate random variables not only using their *parameters*, but also using their *moments*, which often useful in the field of Structural and System Reliability Analysis.

```@setup GenerateRV
using Fortuna
```

## Generating Random Variables Using Moments

To generate a random variable using its moments use pass `"M"` or `"Moments"` as the second argument of `generaterv()` function followed by the moments themselves. 

```@example GenerateRV
# Generate a lognormally distributed random variable R 
# with mean (μ) of 15 and standard deviation (σ) of 10:
R = generaterv("LogNormal", "M", [15, 10])
println("μ = $(mean(R))")
println("σ = $(std(R))")
```

## Generating Random Variables Using Parameters

To generate a random variable using its parameters use pass `"P"` or `"Parameters"` as the second argument of `generaterv()` function followed by the parameters themselves. 

```@example GenerateRV
# Generate a gamma-distributed random variable Q 
# with shape parameter (α) of 16 and scale parameter (θ) of 0.625:
Q = generaterv("Gamma", "P", [16, 0.625])
println("α = $(params(Q)[1])")
println("θ = $(params(Q)[2])")
```

## Supported Random Variables

!!! note
    If you want to define a random variable that is not supported by Fortuna.jl package, please raise an issue on the [Github Issues](https://github.com/AkchurinDA/Fortuna.jl/issues) page.

`Fortuna.jl` package currently supports the following distributions:
- [Exponential](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Exponential)
- [Gamma](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Gamma)
- [Gumbel](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Gumbel)
- [LogNormal](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.LogNormal)
- [Normal](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Normal)
- [Poisson](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Poisson)
- [Uniform](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Uniform)
- [Weibull](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Weibull)

### API

```@docs
generaterv(DistributionName::String, DefineBy::String, Values::Union{Real, Vector{<:Real}})
```
