
module QLawIndirectOptimization

#using Suppressor
using StaticArrays
using AstroUtils, AstroEOMs
using DelimitedFiles
using DifferentialEquations
#using DiffEqSensitivity
#using DrWatson
using LinearAlgebra
#using Zygote
#using NLopt
#using NLsolve
#using ForwardDiff
#using DiffResults
using Interpolations
import MAT
#import Convex, SCS, COSMO
#import Heuristics
import GlobalOptimization
#const MOI = Convex.MOI

# Plotting
#import GLMakie as GM
import CairoMakie as CM

# Define some functions for easy use of 
# c generated code
pow(a,b)    = a^b
fabs(a)     = abs(a)
fmod(a,b)   = mod(a,b)

# Include qlaw source code
include("./qLaw/qLawParams.jl")
#include("./qLaw/qLawThrustAngles.jl")
include("./qLaw/dQn.jl")
include("./qLaw/Qpartials_keplerian.jl")
include("./qLaw/gaussVarKeplerian.jl")
#include("./qLaw/qLawEffectivity_Keplerian.jl")
#include("./qLaw/qLawCoastContinuousCallbackCheck.jl")
#include("./qLaw/qLawSolvers/quickestDescentSolve.jl")
#include("./qLaw/qLawSolvers/quickestDescentSunAngleAnalytic.jl")
#include("./qLaw/qLawSolvers/steepestDescentControl.jl")
include("./qLaw/qLawEOMsSundmanTransformedZOH.jl")
#include("./qLaw/qLawOriginal.jl")
include("./qLaw/qLawThrust_Keplerian.jl")
include("./qLaw/qLaw_transfer_cache.jl")
include("./qLaw/qLaw_transfer_generation.jl")

# Include indirect optimization source code
# include("./indirectOptimization/memfNonlinearFunction.jl")
# include("./indirectOptimization/memfSolve.jl")
# include("./indirectOptimization/minFuelMayerNonlinearFunction.jl")
# include("./indirectOptimization/minFuelMayerSolve.jl")
# include("./indirectOptimization/FinalStateInterpolant.jl")
# include("./indirectOptimization/minFuelMayerPseudoarclengthNonlinearFunction.jl")
# include("./indirectOptimization/pcContinuation.jl")
# include("./indirectOptimization/ACT/rotationMatrix.jl")
# include("./indirectOptimization/ACT/adjointControlTransformation.jl")

export qLawParams, qLawThrustAngles, dQn, qLawCoastContinuousCallbackCheck
export qLawOriginal
export qLawThrust_Keplerian
export generate_qlaw_transfer
export plot_transfer, plot_transfer_gl
export dump_to_mat
#export memfSolve, minFuelMayerSolve
#export adjointControlTransformation
#export FinalStateInterpolant, getFinalStateAndTime
#export minFuelMayerPseudoarclengthNonlinearFunction!
#export pcContinuation

end