# Sensitivity Analysis

## Overview

```math
\nabla_{\vec{\theta}_{g}} \beta = \dfrac{\nabla_{\vec{\theta}_{g}} g(\vec{x}^{*}, \vec{\theta}_{g})}{||\nabla_{\vec{u}} G(\vec{u}^{*}, \vec{\theta}_{g})||}
```

```math
\nabla_{\vec{\theta}_{g}} P_{f} = -\phi(\beta) \nabla_{\vec{\theta}_{g}} \beta
```

## API

```@docs
SensitivityProblem
analyze(Problem::SensitivityProblem)
```