
module QLawIndirectOptimization

using StaticArrays
using AstroUtils, AstroEOMs
using DelimitedFiles
using DifferentialEquations
using LinearAlgebra
using Interpolations
using LaTeXStrings

import MAT
import GlobalOptimization
import CairoMakie as CM

# Define some functions for easy use of 
# c generated code
pow(a,b)    = a^b
fabs(a)     = abs(a)
fmod(a,b)   = mod(a,b)

# Include qlaw source code
include("qLawParams.jl")
include("utils.jl")
include("Qpartials_keplerian.jl")
include("gaussVarKeplerian.jl")
include("qLawEOMsSundmanTransformedZOH.jl")
include("qLawThrust_Keplerian.jl")
include("qLaw_transfer_cache.jl")
include("qLaw_transfer_generation.jl")
include("publication_plotting.jl")

export qLawParams
export generate_qlaw_transfer
export plot_transfer
export dump_to_mat
export generate_publication_figures

end