module STRELA
# Load dependencies:
using Distributions: Distribution
using Distributions: MvNormal
using Distributions: Exponential, Gamma, Gumbel, LogNormal, Normal, Poisson, Uniform
using Distributions: mean, std, pdf, cdf, quantile
using FastGaussQuadrature: gausslegendre
using ForwardDiff: gradient
using LinearAlgebra
using Random: rand, randn, shuffle
using NonlinearSolve: NonlinearProblem
using NonlinearSolve: NewtonRaphson
using NonlinearSolve: solve

# Include the following files into the scope of the module:
include("DefineStructures.jl")
export NatafTransformation
export ReliabilityProblem
include("GenerateRandomVariables.jl")
export generaterv
export convertmoments
include("SampleRandomVariables.jl")
export samplerv
include("PerformNatafTransformation.jl")
export getdistortedcorrelation
export transformsamples
export getjacobian
export jointpdf
include("PerformReliabilityAnalysis.jl")
export MCFOSM
export FORM
end
