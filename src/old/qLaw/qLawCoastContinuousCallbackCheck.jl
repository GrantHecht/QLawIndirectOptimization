
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
    ηa      = dQ / dQmin

    val     = 1.0
    if ps.ηr > 0.0
        val = ηr - ps.ηr
    elseif ps.ηa > 0.0
        val = ηa - ps.ηa
    end

    return val
end