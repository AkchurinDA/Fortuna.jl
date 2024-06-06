module Fortuna
# --------------------------------------------------
# IMPORT PACKAGES
# --------------------------------------------------
import Base
import Distributions
import FastGaussQuadrature
import FiniteDiff
import ForwardDiff
import LinearAlgebra
import NonlinearSolve
import Random
import SpecialFunctions
using  DocStringExtensions

import ProgressMeter

# --------------------------------------------------
# REEXPORT PACKAGES
# --------------------------------------------------
import Reexport
Reexport.@reexport import Distributions: rand, pdf              # Extended functions
Reexport.@reexport import Distributions: mean, std, cor, params # Useful functions
Reexport.@reexport import LinearAlgebra: I

# --------------------------------------------------
# DEFINE ABSTRACT TYPES
# --------------------------------------------------
"""
    abstract type AbstractIsoprobabilisticTransformation end

Abstract type for isoprobabilistic transformations.
"""
abstract type AbstractIsoprobabilisticTransformation end

# """
#     abstract type AbstractSamplingTechnique end

# Abstract type for sampling techniques.
# """
# abstract type AbstractSamplingTechnique end

"""
    abstract type AbstractReliabilityProblem end

Abstract type for reliability problems.
"""
abstract type AbstractReliabilityProblem end

"""
    abstract type AbstractReliabililyAnalysisMethod end

Abstract type for reliability analysis methods.
"""
abstract type AbstractReliabililyAnalysisMethod end

"""
    abstract type FORMSubmethod end

Abstract type for First-Order Reliability Method's (FORM) submethods.
"""
abstract type FORMSubmethod end

"""
    abstract type SORMSubmethod end

Abstract type for Second-Order Reliability Method's (FORM) submethods.
"""
abstract type SORMSubmethod end

# --------------------------------------------------
# EXPORT TYPES AND FUNCTIONS
# --------------------------------------------------
include("Isoprobabilistic Transformations/NatafTransformation.jl")
include("Isoprobabilistic Transformations/RosenblattTransformation.jl")
include("Random Variables/DefineRandomVariables.jl")
include("Random Variables/SampleRandomVariables.jl")
include("Reliability Problems/ReliabilityProblems.jl")
include("InverseReliabilityProblems.jl")
include("SensitivityProblems.jl")
# include("Utilities/MakieRecipes.jl") # Temporarily removed to avoid dependency on Makie - need to create an extension.
export AbstractSamplingTechnique
export ITS, LHS
export AbstractTransformation
export NatafTransformation, RosenblattTransformation
export AbstractReliabilityProblem
export ReliabilityProblem
export MC, MCCache
export IS, ISCache
export FORM, MCFOSM, MCFOSMCache, RF, RFCache, HLRF, HLRFCache, iHLRF, iHLRFCache
export SORM, CF, CFCache, PF, PFCache
export SSM, SSMCache
export SensitivityProblemTypeI, SensitivityProblemTypeII, SensitivityProblemCache
export InverseReliabilityProblem, InverseReliabilityProblemCache
export randomvariable
export getdistortedcorrelation, transformsamples, getjacobian
export solve
end