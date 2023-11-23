# Second-Order Reliability Methods

The SORM is an improvement over the FORM by accounting for the curved nature of the failure boundary given by ``g(\vec{X}) = 0`` around the design point ``\vec{x}^{*}``; thus, providing a better approximation of the probability of failure ``P_{f}``.

## Curve-Fitting Method

`Fortuna.jl` package implements the CF method that fits a hyper-paraboloid surface with a vertex at the design point ``\vec{x}^{*}`` and the principal curvatures matching the principal curvatures of the failure boundary  given by ``g(\vec{X}) = 0`` at that point. The probabilities ``P_{f}`` of failure are estimated using [Hohenbichler1988](@citet) and [Breitung1984](@citet) approximations of the exact solution provided by [Tvedt1990](@citet). The calculated probabilities of failure ``P_{f}`` are then used to estimate the generalized reliability indices ``\beta``, which account for the curved nature of the failure boundary given by ``g(\vec{X}) = 0`` around the design point ``\vec{x}^{*}``.

## Gradient-Free Method

!!! note
    This feature is currently under development.


## Point-Fitting Method

!!! note
    This feature is currently under development.
