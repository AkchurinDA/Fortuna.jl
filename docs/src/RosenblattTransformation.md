# Rosenblatt Transformation

!!! note
    This feature is currently under development.

## Overview

The Rosenblatt Transformation is another widely utilized isoprobabilistic transformation in structural reliability analysis. Similar to the Nataf Transformation, its purpose is to transform random vectors with correlated non-normal marginal random variables ``\vec{X}`` into random vectors with uncorrelated standard normal marginal random variables ``\vec{U}``. Murray Rosenblatt introduced this transformation in 1952 [Rosenblatt1952](@cite).

The Nataf Transformation ``\vec{U} = T_{R}(\vec{X})`` is composed of two transformations: 

```math
\vec{U} = T_{N}(\vec{X}) = (T_{1} \circ T_{2})(\vec{X})
```

- The first transformation ``\vec{Z} = T_{1}(\vec{X})`` transforms random vector with *correlated non-normal marginal random variables* ``\vec{X}`` (with correlation matrix ``\rho^{X}``) into random vector with *uncorrelated uniform marginal random variables* ``\vec{Z}``.

```math
\vec{Z} = T_{1}(\vec{X}) = \begin{bmatrix} F_{X_{1}}(X_{1}) \\ \vdots \\ F_{X_{n} | X_{n - 1}, \dots, X_{1}}(X_{n} | X_{n - 1}, \dots, X_{1}) \end{bmatrix}
```

- The second transformation ``\vec{U} = T_{2}(\vec{Z})`` transforms random vector with *uncorrelated uniform marginal random variables* ``\vec{Z}`` into random vector with *uncorrelated standard normal marginal random variables* ``\vec{U}``.

```math
\vec{U} = T_{2}(\vec{Z}) = \begin{bmatrix} \Phi^{-1}(Z_{1}) \\ \vdots \\ \Phi^{-1}(Z_{n}) \end{bmatrix}
```

## Associated Types and Functions

```@docs
RosenblattTransformation
```