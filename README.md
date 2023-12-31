<div align="center">
  <img src="assets/logo.svg" alt = "Logo" width="35%">

  | Contributors | [Damir Akchurin](https://scholar.google.com/citations?user=chYaDcIAAAAJ&hl=en) |
  | :---: | :---: |
  | Build Status | [![Build Status](https://github.com/AkchurinDA/Fortuna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AkchurinDA/Fortuna.jl/actions/workflows/CI.yml) |
  | Documentation | [![Documentation](https://img.shields.io/badge/Documentation-Stable-blue.svg)](https://AkchurinDA.github.io/Fortuna.jl/stable) |
  | Latest Release | [![Latest Release](https://juliahub.com/docs/Fortuna/version.svg)](https://github.com/AkchurinDA/Fortuna.jl/releases) |
  | Downloads | [![Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Fortuna&label=Downloads)](https://pkgs.genieframework.com?packages=Fortuna) |
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

For assistance with the package, please raise an issue on the Github Issues page. Please use the appropriate labels to indicate the specific functionality you are inquiring about. Alternatively, contact the author directly.

## Acknowledgements

The author thanks the academic and industrial partners of the Cold-Formed Steel Research Consortium’s [“Reliability 2030”](https://cfsrc.org/2023/01/01/reliability-2030-design-of-steel-as-a-system/) initiative for their financial support.

## Roadmap

The following functionality is planned to be added:
- [x] Isoprobabilistic Transformations
  - [x] Nataf Transformation
  - [ ] Rosennblatt Transformation
- [x] Sampling Techniques
  - [x] Inverse Transform Sampling
  - [x] Latin Hypercube Sampling
- [ ] Reliability Methods
  - [ ] First-Order Reliability Method
    - [x] Mean-Centered First-Order Second-Moment Method
    - [ ] Hasofer-Lind Method
    - [ ] Rackwitz-Fiessler Method
    - [x] Hasofer-Lind Rackwitz-Fiessler Method
    - [x] Improved Hasofer-Lind Rackwitz-Fiessler Method
  - [ ] Second-Order Reliability Method
    - [x] Curve-Fitting Method
    - [ ] Gradient-Free Method
    - [ ] Point-Fitting Method
  - [x] Subset Simulation Method
- [ ] Sensitivity Analysis