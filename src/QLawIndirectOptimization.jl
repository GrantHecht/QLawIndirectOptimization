
module QLawIndirectOptimization 
    using StaticArrays

    # Define some functions for easy use of 
    # c generated code
    pow(a,b)    = a^b
    fabs(a)     = abs(a)
    fmod(a,b)   = mod(a,b)

    include("qLawParams.jl")
    include("qLawThrust.jl")
    include("dQn.jl")
    include("QLawEOMs.jl")
    include("qLaw.jl")

    export qLawParams, qLaw, qLawCoastContinuousCallbackCheck
    export qLawThrustAngles, dQn
end