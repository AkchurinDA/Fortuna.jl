module Fortuna
# Reexport Distributions package:
using Reexport
@reexport using Distributions

# Load dependencies:
using DocStringExtensions
using FastGaussQuadrature: gausslegendre
using ForwardDiff: gradient, hessian
using LinearAlgebra
using Random: rand, randn, shuffle
using NonlinearSolve: NonlinearProblem, IntervalNonlinearProblem
using NonlinearSolve: NewtonRaphson, Bisection
using NonlinearSolve: solve
using SpecialFunctions: gamma

# Include the following files into the scope of the module:
include("Types.jl")
export AbstractTransformation
export NatafTransformation, RosenblattTransformation
export AbstractSamplingTechnique
export ITS, LHS
export AbstractReliabilityProblem
export ReliabilityProblem
export AbstractReliabililyAnalysisMethod
export FORMSubmethod, SORMSubmethod
export MCS, MCSCache
export FORM, MCFOSM, HLRF, iHLRF, MCFOSMCache, HLRFCache, iHLRFCache
export SORM, CF, PF, CFCache, PFCache
export SSM, SSMCache
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
