# News

## Release V0.5.1

- Added measure of importance of random variables $\gamma$ to FORM. Note that the importance vector $\gamma$ is only available if the reliability analysis was carried out using HLRF or iHLRF methods.

```julia
# Plain Hasofer-Lind-Rackwitz-Fiessler method:
FORMSolution = analyze(Problem, FORM(HLRF()))
FORMSolution.γ

# Improved Hasofer-Lind-Rackwitz-Fiessler method:
FORMSolution = analyze(Problem, FORM(iHLRF()))
FORMSolution.γ
```

## Release V0.5.0

- Added functionality to probabilities of failure by means of Monte Carlo Simulations.

```julia
MCSSolution = analyze(Problem, MCS())
```

- Added Point-Fitting method that falls within a broader class of the Second-Order Reliability Methods (SORM).

```julia
SORMSolution = analyze(Problem, SORM(PF()))
```

- Added an option to sample random variables with correlated marginal random variables $\vec{X}$ using any implemented sampling technique. Note that your choice of the sampling technique pertains to the sampling technique that is used to generate samples in the space of uncorrelated standard normal random variables $\vec{U}$. The generated samples are then transformed from $\vec{U}$- to $\vec{X}$-space according to transformation you use to define a random vector with correlated marginal random variables $\vec{X}$, i.e., either Nataf or Rosennblatt (not yet implemented) transformation.

```julia
XSamples, ZSamples, USamples = samplerv(TransformationObject, 1000, ITS())
XSamples, ZSamples, USamples = samplerv(TransformationObject, 1000, LHS())
```

- Updated the documentation.
- Added Dependabot to automatically take care of GitHub Actions' versions.

## Release V0.4.1

- Accelerate the Subset Simulation Method.
- 
```julia
SSMSolution = analyze(Problem, SSM())
```

## Release V0.4.0

- Added the Subset Simulation Method to compute probabilities of rare failure events.

## Release V0.3.4

- General housekeeping and bug hunting.
- Added a simple documentation describing the basic functionality of the package.
- `Fortuna.jl` package now reexports `Distributions.jl` package. This allows to access basic functions, like `mean()`, `std()`, `params()`, etc., without loading `Distributions.jl` package separately.

## Release V0.3.3

- Removed `println()` statement from `samplerv()` function that was left there by accident.

## Release V0.3.2

- Fixed the generalized reliability index calculation in the SORM.

## Release V0.3.1

- Added error-catching for the SORM's probability of failure approximations.

## Release V0.3.0

- Added Curve-Fitting method that falls within a broader class of the Second-Order Reliability Methods (SORM).
- Added cache output system for analysis results for easier work in the future. You can now access a lot of iterative data related to the reliability analysis, such as the value of the limit state function, its gradient, merit function, and more at each iteration.

```julia
FORMSolution = analyze(Problem, FORM())
SORMSolution = analyze(Problem, SORM())
```

- Benchmarked FORM, SORM, and Nataf Transformation functionalities against various sources. The sources for the benchmarks are indicated. 
- Added support for Weibull distribution.

## Release V0.2.0

- Sampling function `samplerv()` now uses native `rand()` and `randn()` functions from the `Random` package to perform ITS and LHS sampling.
- `SamplingTechnique` argument's options in `samplerv()` function have been changed from strings (`"ITS"`, `"LHS"`) to their own types (`ITS()`, `LHS()`).
- Added `analyze()` function that replaces `MCFOSM()` and `FORM()` functions. The deprecated functions `MCFOSM()` and `FORM()` are now `AnalysisMethod` argument's options in `analyze()` function with their own types. For example, to perform reliability analysis using FORM use the following syntax:

```julia
β       = analyze(Problem, FORM(MCFOSM()))
β, x, u = analyze(Problem, FORM(HLRF()))
β, x, u = analyze(Problem, FORM(iHLRF()))
```

## Release V0.1.0

- Initial release of `Fortuna.jl`.