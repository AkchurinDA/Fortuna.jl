# Nataf Transformation

## Theory

```math
\underline{U} = T(\underline{X}) = (T_{1} \circ T_{2})(\underline{X})
```

```math
\underline{Z} = T_{1}(\underline{X}) = \begin{bmatrix} \Phi^{-1}(F_{X_{1}}(x_{1})) \\ \vdots \\ \Phi^{-1}(F_{X_{n}}(x_{n})) \end{bmatrix}
```

```math
\underline{U} = T_{2}(\underline{Z}) = \underline{\underline{\Gamma}} \ \underline{Z}
```

## Associated Functions

```@docs
getdistortedcorrelation(X::Vector{<:Distribution}, ρˣ::Matrix{<:Real})
jointpdf(Object::NatafTransformation, XSamples::Union{Vector{<:Real},Matrix{<:Real}})
```