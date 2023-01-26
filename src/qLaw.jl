
# Include qlaw source code
include(srcdir("qLawParams.jl"))
include(srcdir("qLawThrustUnitVector.jl"))
include(srcdir("dQn.jl"))

function qLaw(oe,ps::qLawParams)
    # Compute qlaw thrust direction
    u = qLawThrustUnitVector(oe[1],oe[2],oe[3],oe[4],oe[5],oe[6],ps)

    # Compute effectivity
    dQmin   = Inf
    dQmax   = -Inf
    dQ      = dQn(oe[1],oe[2],oe[3],oe[4],oe[5],oe[6],ps)
    θs      = range(0.0,2*pi; length = ps.steps)
    for i in eachindex(θs)
        dQθ = dQn(oe[1],oe[2],oe[3],oe[4],oe[5],θs[i],ps)
        if dQθ < dQmin; dQmin = dQθ; end
        if dQθ > dQmax; dQmax = dQθ; end
    end
    ηr      = (dQ - dQmax) / (dQmin - dQmax)
    coast   = ηr < ps.ηr

    return (u,coast)
end

function qLawCoastFlagCheck(oe,ps)
     # Compute qlaw thrust direction
    u = qLawThrustUnitVector(oe[1],oe[2],oe[3],oe[4],oe[5],oe[6],ps)

    # Compute effectivity
    dQmin   = Inf
    dQmax   = -Inf
    dQ      = dQn(oe[1],oe[2],oe[3],oe[4],oe[5],oe[6],ps)
    θs      = range(0.0,2*pi; length = ps.steps)
    for i in eachindex(θs)
        dQθ = dQn(oe[1],oe[2],oe[3],oe[4],oe[5],θs[i],ps)
        if dQθ < dQmin; dQmin = dQθ; end
        if dQθ > dQmax; dQmax = dQθ; end
    end
    ηr      = (dQ - dQmax) / (dQmin - dQmax)
    coast   = ηr < ps.ηr

    swapFlag = false
    if coast != ps.coasting
        swapFlag = true
    end
    return swapFlag
end

function qLawCoastContinuousCallbackCheck(oe,ps)
     # Compute qlaw thrust direction
    u = qLawThrustUnitVector(oe[1],oe[2],oe[3],oe[4],oe[5],oe[6],ps)

    # Compute effectivity
    dQmin   = Inf
    dQmax   = -Inf
    dQ      = dQn(oe[1],oe[2],oe[3],oe[4],oe[5],oe[6],ps)
    θs      = range(0.0,2*pi; length = ps.steps)
    for i in eachindex(θs)
        dQθ = dQn(oe[1],oe[2],oe[3],oe[4],oe[5],θs[i],ps)
        if dQθ < dQmin; dQmin = dQθ; end
        if dQθ > dQmax; dQmax = dQθ; end
    end
    ηr      = (dQ - dQmax) / (dQmin - dQmax)

    return ηr - ps.ηr
end