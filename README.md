<div align="center">
  <img src="assets/social-preview.svg" alt = "Social Preview" width = 75%>

  **Audentes Fortuna Iuvat**

  | Developer | [Damir Akchurin](https://scholar.google.com/citations?user=chYaDcIAAAAJ&hl=en) |
  | :--- | :--- |
  | Build Status | [![Build Status](https://github.com/AkchurinDA/Fortuna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AkchurinDA/Fortuna.jl/actions/workflows/CI.yml) |
  | Latest Release | [![Latest Release](https://juliahub.com/docs/Fortuna/version.svg)](https://github.com/AkchurinDA/Fortuna.jl/releases) |
  | Documentation | [![Documentation](https://img.shields.io/badge/Documentation-Stable-blue.svg)](https://AkchurinDA.github.io/Fortuna.jl/stable) <br> [![Documentation](https://img.shields.io/badge/Documentation-Dev-blue.svg)](https://AkchurinDA.github.io/Fortuna.jl/dev) |
  | Downloads | [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FFortuna&query=total_requests&label=Total)](http://juliapkgstats.com/pkg/Fortuna) <br> [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FFortuna&query=total_requests&label=Monthly&suffix=%2FMonth)](http://juliapkgstats.com/pkg/Fortuna) |
  | License | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/AkchurinDA/Fortuna.jl/blob/main/LICENSE.md) |

</div>

## Description

`Fortuna.jl` is a general-purpose Julia package for structural and system reliability analysis.

## Installation

To install `Fortuna.jl` package, type `]` in Julia REPL to enter the built-in Julia package manager and execute the following command:

```
pkg> add Fortuna
```

## License

`Fortuna.jl` package is distributed under the [MIT license](https://en.wikipedia.org/wiki/MIT_License). More information can be found in the [`LICENSE.md`](https://github.com/AkchurinDA/Fortuna.jl/blob/main/LICENSE.md) file.

## Help and Support

For assistance with the package, please raise an issue on the [Github Issues](https://github.com/AkchurinDA/Fortuna.jl/issues) page. Please use the appropriate labels to indicate the specific functionality you are inquiring about. Alternatively, contact the author directly at [AkchurinDA@gmail.com](mailto:AkchurinDA@gmail.com?subject=Fortuna.jl).

## Acknowledgements

The author thanks the academic and industrial partners of the [“Reliability 2030”](https://cfsrc.org/2023/01/01/reliability-2030-design-of-steel-as-a-system/) initiative for their financial support.

## Roadmap

The following functionality is planned to be added:

- [x] Sampling techniques
    - [x] Inverse transform sampling
    - [x] Latin hypercube sampling
- [x] Isoprobabilistic transformations
    - [x] Nataf transformation
- [x] Reliability analysis
  - [x] Monte Carlo methods
    - [x] Direct Monte Carlo simulations
    - [x] Importance sampling method
  - [x] First-order reliability methods
    - [x] Mean-centered first-order second-moment method
    - [x] Rackwitz-Fiessler method
    - [x] Hasofer-Lind Rackwitz-Fiessler method
    - [x] Improved Hasofer-Lind Rackwitz-Fiessler method
  - [x] Second-order reliability methods
    - [x] Curve-fitting method
    - [x] Point-fitting method
  - [x] Subset simulation method
- [x] Inverse reliability analysis
- [x] Sensitivity analysis
  - [x] w.r.t. parameters of limit state functions
  - [x] w.r.t. parameters and moments of random vectors
- [x] Limit state functions defined using
  - [x] "Internal" FE models
  - [x] "External" FE models
  - [x] Surrogate models