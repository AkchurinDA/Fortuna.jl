# Rosenblatt Transformation

## Theory

```math
\begin{align*}
\underline{U} &= T_{R}(\underline{X}) = (T_{1} \circ T_{2})(\underline{X}) \\
\underline{Z} &= T_{1}(\underline{X}) = \begin{bmatrix} F_{X_{1}}(X_{1}) \\ \vdots \\ F_{X_{n} | X_{1}, \dots, X_{n - 1}}(X_{n} | X_{1}, \dots, X_{n - 1}) \end{bmatrix} \\
\underline{U} &= T_{2}(\underline{Z}) = \begin{bmatrix} \Phi^{-1}(Z_{1}) \\ \vdots \\ \Phi^{-1}(Z_{n}) \end{bmatrix}
\end{align*}
```