
function qLawCoastContinuousCallbackCheck(oe,m,ps)
    # Compute effectivity
    dQmin   = Inf
    dQmax   = -Inf
    dQ      = dQn(oe[1],oe[2],oe[3],oe[5],oe[4],oe[6],m,ps)
    θs      = range(0.0,2*pi; length = ps.steps)
    for i in eachindex(θs)
        dQθ = dQn(oe[1],oe[2],oe[3],oe[5],oe[4],θs[i],m,ps)
        if dQθ < dQmin; dQmin = dQθ; end
        if dQθ > dQmax; dQmax = dQθ; end
    end
    ηr      = (dQ - dQmax) / (dQmin - dQmax)

    return ηr - ps.ηr
end