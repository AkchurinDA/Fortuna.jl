# First-Order Reliability Methods

## Mean-Centered First-Order Second-Moment Method

The MCFOSM method is the simplest and least computationally expensive type of reliability method. It utilizes the first-order Taylor expansion of the limit state function ``g(\vec{X})`` at the mean values and the first two moments of the marginal random variables involved in the reliability problem to evaluate the reliability index. However, despite the fact that it is simple and does not require the complete knowledge of the random variables involved in the reliability problem, the MCFOSM method faces an issue known as the *invariance problem*. This problem arises because the resulting reliability index ``\beta`` is dependent on the formulation of the limit state function ``g(\vec{X})``. In other words, two equivalent limit state functions with the same failure boundaries produce two different reliability indices; thus, the use of MCFOSM method is not recommended.

## Hasofer-Lind Method

!!! note
    This feature is currently under development.

## Rackwitz-Fiessler Method

!!! note
    This feature is currently under development.

## Plain and Improved Hasofer-Lind-Rackwitz-Fiessler Method

The HLRF method overcomes the invariance problem faced by the MCFOSM method by using the first-order Taylor expansion of the limit state function at a point known as the *design point* ``\vec{x}^{*}`` that lies on the failure boundary given by ``g(\vec{X}) = 0``. Since the design point is not known a priori, the HLRF method is inherently an iterative method. `Fortuna.jl` implements two versions of HLRF method: plain HLRF method where the step size in the negative gradient descent is always set to unity and improved HLRF (iHLRF) method where the step size is determined using a line search algorithm.
