module Fortuna
# Load dependencies:
using Distributions: Distribution
using Distributions: MvNormal
using Distributions: Exponential, Gamma, Gumbel, LogNormal, Normal, Poisson, Uniform
using Distributions: mean, std, pdf, cdf, quantile, cor
using FastGaussQuadrature: gausslegendre
using ForwardDiff: gradient
using LinearAlgebra
using Random: rand, randn, shuffle
using NonlinearSolve: NonlinearProblem
using NonlinearSolve: NewtonRaphson
using NonlinearSolve: solve

# Include the following files into the scope of the module:
include("Structures.jl")
export NatafTransformation, RosenblattTransformation
export ITS, LHS
export ReliabilityProblem
export MCFOSM
export FORM, HLRF, iHLRF
include("GenerateRandomVariables.jl")
export generaterv
export convertmoments
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
