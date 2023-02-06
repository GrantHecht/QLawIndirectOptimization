
module QLawIndirectOptimization

using StaticArrays
using AstroUtils, AstroEOMs
using DelimitedFiles
using DifferentialEquations
using DrWatson
using LinearAlgebra
using Zygote
using NLopt
import Convex, SCS
const MOI = Convex.MOI

# Define some functions for easy use of 
# c generated code
pow(a,b)    = a^b
fabs(a)     = abs(a)
fmod(a,b)   = mod(a,b)

# Include qlaw source code
include("./qLaw/qLawParams.jl")
include("./qLaw/qLawThrustAngles.jl")
include("./qLaw/dQn.jl")
include("./qLaw/Qpartials_keplerian.jl")
include("./qLaw/qLawCoastContinuousCallbackCheck.jl")
include("./qLaw/qLawSolvers/quickestDescentSolve.jl")
include("./qLaw/qLawSolvers/steepestDescentControl.jl")
include("./qLaw/qLawEOMsSundmanTransformedZOH.jl")
include("./qLaw/qLawOriginal.jl")
include("./qLaw/qLawThrust_Keplerian.jl")

export qLawParams, qLawThrustAngles, dQn, qLawCoastContinuousCallbackCheck
export qLawOriginal
export qLawThrust_Keplerian

end