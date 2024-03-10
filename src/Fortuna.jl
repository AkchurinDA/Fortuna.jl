module Fortuna
import Base
import Random
import Distributions
import FastGaussQuadrature
import LinearAlgebra
import NonlinearSolve
import SpecialFunctions
import ForwardDiff

using Reexport
# Extended functions:
@reexport import Distributions: rand, pdf 
# Useful functions:
@reexport import Distributions: mean, std
@reexport import LinearAlgebra: I

include("Types.jl")
export AbstractSamplingTechnique
export ITS, LHS
export AbstractTransformation
export NatafTransformation, RosenblattTransformation
export AbstractReliabilityProblem
export ReliabilityProblem
export MC, MCCache
export IS, ISCache
export FORM, MCFOSM, MCFOSMCache, HL, HLCache, RF, RFCache, HLRF, HLRFCache, iHLRF, iHLRFCache
export SORM, CF, CFCache, PF, PFCache
export SSM, SSMCache
export SensitivityProblem, SensitivityProblemCache
export InverseReliabilityProblem, InverseReliabilityProblemCache
include("Random Variables/GenerateRandomVariables.jl")
include("Random Variables/SampleRandomVariables.jl")
export randomvariable
include("Isoprobabilistic Transformations/NatafTransformation.jl")
export getdistortedcorrelation, transformsamples, getjacobian
include("Reliability Problems/ReliabilityProblems.jl")
include("Inverse Reliability Problems/InverseReliabilityProblems.jl")
include("Sensitivity Problems/SensitivityProblems.jl")
export solve
end