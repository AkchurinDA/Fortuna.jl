module Fortuna
# Load dependencies:
using Distributions: Distribution
using Distributions: MvNormal
using Distributions: Exponential, Gamma, Gumbel, LogNormal, Normal, Poisson, Uniform, Weibull
using Distributions: mean, std, pdf, cdf, quantile, cor
using FastGaussQuadrature: gausslegendre
using ForwardDiff: gradient, hessian
using LinearAlgebra
using Random: rand, randn, shuffle
using NonlinearSolve: NonlinearProblem, IntervalNonlinearProblem
using NonlinearSolve: NewtonRaphson, Bisection
using NonlinearSolve: solve
using SpecialFunctions: gamma
using QuadGK: quadgk

# Include the following files into the scope of the module:
include("Structures.jl")
export NatafTransformation, RosenblattTransformation
export ITS, LHS
export ReliabilityProblem
export FORM, MCFOSM, HLRF, iHLRF
export SORM, CF, PF
include("GenerateRandomVariables.jl")
export generaterv
include("SampleRandomVariables.jl")
export samplerv
include("Transformations/NatafTransformation.jl")
export getdistortedcorrelation
export transformsamples
export getjacobian
export jointpdf
include("PerformReliabilityAnalysis.jl")
export analyze
end
