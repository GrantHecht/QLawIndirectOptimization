# Must have project active before including
using StaticArrays

# Include symbolics generated functions (not using this functionallity)
#qLawThrustAngles = include(srcdir("qLawThrustAngles.jl"))
#Qn = include(srcdir("Qn.jl"))

# Include source
include(srcdir("qLaw.jl"))