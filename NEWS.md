# News

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