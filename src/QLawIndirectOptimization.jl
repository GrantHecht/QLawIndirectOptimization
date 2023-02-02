
module QLawIndirectOptimization

using StaticArrays
using AstroUtils, AstroEOMs
using DelimitedFiles
using DifferentialEquations
using DrWatson

# Define some functions for easy use of 
# c generated code
pow(a,b)    = a^b
fabs(a)     = abs(a)
fmod(a,b)   = mod(a,b)

# Include qlaw source code
include("qLawParams.jl")
include("qLawThrustAngles.jl")
include("dQn.jl")
include("qLawCoastContinuousCallbackCheck.jl")
include("qLawEOMsSundmanTransformedZOH.jl")
include("qLaw.jl")

export qLawParams, qLawThrustAngles, dQn, qLawCoastContinuousCallbackCheck
export qLaw

end