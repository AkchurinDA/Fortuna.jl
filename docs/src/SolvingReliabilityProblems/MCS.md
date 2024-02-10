# Monte Carlo Simulations

## Overview

```math
P_{f} = P(\Omega_{f}) = \int_{\Omega_{f}} f_{\vec{X}}(\vec{x}) d\vec{x} = \int_{\mathbb{R}^{n}} \mathbb{I}(\vec{x}) f_{\vec{X}}(\vec{x}) d\vec{x} = \mathbb{E}[\mathbb{I}(\vec{x})]
```

```math
\mathbb{I}(\vec{x}) = 
\begin{cases}
    1 & \text{if } \vec{x} \in \Omega_{f} \\
    0 & \text{otherwise}
\end{cases}
```

```math
P_{f}^{MCS} = \dfrac{1}{N} \sum_{i = 1}^{N} \mathbb{I}(\vec{x}_{i})
```

### API