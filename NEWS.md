# News

# Release V0.4.1
- Accelerate the Subset Simulation Method.

# Release V0.4.0
- Added the Subset Simulation Method to compute probabilities of rare failure events.

# Release V0.3.4
- General housekeeping and bug hunting.
- Added a simple documentation describing the basic functionality of the package.
- `Fortuna.jl` package now reexports `Distributions.jl` package. This allows to access basic functions, like `mean()`, `std()`, `params()`, etc., without loading `Distributions.jl` package separately.

# Release V0.3.3
- Removed `println()` statement from `samplerv()` function that was left there by accident.

# Release V0.3.2
- Fixed the generalized reliability index calculation in the SORM.

# Release V0.3.1
- Added error-catching for the SORM's probability of failure approximations.

# Release V0.3.0
- Added Curve-Fitting method that falls within a broader class of the Second-Order Reliability Methods (SORM).
- Added cache output system for analysis results for easier work in the future. You can now access a lot of iterative data related to the reliability analysis, such as the value of the limit state function, its gradient, merit function, and more at each iteration.

```julia
FORMSolution = analyze(Problem, FORM())
SORMSolution = analyze(Problem, SORM())
```

- Benchmarked FORM, SORM, and Nataf Transformation functionalities against various sources. The sources for the benchmarks are indicated. 
- Added support for Weibull distribution.

# Release V0.2.0
- Sampling function `samplerv()` now uses native `rand()` and `randn()` functions from the `Random` package to perform ITS and LHS sampling.
- `SamplingTechnique` argument's options in `samplerv()` function have been changed from strings (`"ITS"`, `"LHS"`) to their own types (`ITS()`, `LHS()`).
- Added `analyze()` function that replaces `MCFOSM()` and `FORM()` functions. The deprecated functions `MCFOSM()` and `FORM()` are now `AnalysisMethod` argument's options in `analyze()` function with their own types. For example, to perform reliability analysis using FORM use the following syntax:

```julia
β = analyze(Problem, FORM(MCFOSM()))
β, x, u = analyze(Problem, FORM(HLRF()))
β, x, u = analyze(Problem, FORM(iHLRF()))
```

# Release V0.1.0
- Initial release of `Fortuna.jl`.