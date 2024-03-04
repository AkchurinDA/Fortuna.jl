module Fortuna
# Reexport some package and their functionalities:
using Reexport
@reexport using Distributions
@reexport using LinearAlgebra: I

# Load dependencies:
using Base.Iterators:       product, repeated
using DocStringExtensions
using FastGaussQuadrature:  gausslegendre
using ForwardDiff:          gradient, hessian
using LinearAlgebra
using Random:               rand, randn, shuffle
using NonlinearSolve:       NonlinearProblem, IntervalNonlinearProblem
using NonlinearSolve:       NewtonRaphson, Bisection
using NonlinearSolve:       solve
using SpecialFunctions:     gamma

# Include the following files into the scope of the module:
include("Types.jl")
export AbstractSamplingTechnique
export ITS, LHS
export AbstractTransformation
export NatafTransformation, RosenblattTransformation
export AbstractReliabilityProblem
export ReliabilityProblem, InverseReliabilityProblem, SensitivityProblem
export SensitivityProblemCache
export AbstractReliabililyAnalysisMethod
export FORMSubmethod, SORMSubmethod
export MCS, MCSCache
export IS, ISCache
export FORM, MCFOSM, HLRF, iHLRF, MCFOSMCache, HLRFCache, iHLRFCache
export SORM, CF, PF, CFCache, PFCache
export SSM, SSMCache
include("Random Variables/GenerateRandomVariables.jl")
include("Random Variables/SampleRandomVariables.jl")
export generaterv
export samplerv
include("Isoprobabilistic Transformations/NatafTransformation.jl")
include("Isoprobabilistic Transformations/RosenblattTransformation.jl")
export getdistortedcorrelation
export transformsamples
export getjacobian
export jointpdf
include("Reliability Problems/ReliabilityProblems.jl")
include("Inverse Reliability Problems/InverseReliabilityProblems.jl")
include("Sensitivity Problems/SensitivityProblems.jl")
export analyze
end
