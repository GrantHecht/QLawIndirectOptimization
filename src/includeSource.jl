# Must have project active before including
using SymbolicUtils, StaticArraysCore

# Include symbolics generated functions
qLawThrustAngles = include(srcdir("qLawThrustAngles.jl"))
Qn = include(srcdir("Qn.jl"))

# Include source
include(srcdir("qLaw.jl"))