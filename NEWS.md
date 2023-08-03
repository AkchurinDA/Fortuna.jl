# News

# Release V0.3.0
- Benchmarked FORM and Nataf Transformation functionalities against [UQPy](https://github.com/SURGroup/UQpy). The sources for the benchmarks are indicated. 
- Added support for Weibull distribution.
- Probability of failure is now one the outputs of the reliability analysis conducted using HLRF or iHLRF methods, such that:
```julia
β, PoF, x, u = analyze(Problem, FORM(HLRF()))
β, PoF, x, u = analyze(Problem, FORM(iHLRF()))
```
- Added Curve-Fitting and Point-Fitting methods that fall within a broader class of Second-Order Reliability Methods (SORM).

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