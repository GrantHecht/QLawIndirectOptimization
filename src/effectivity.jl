
function effectivity(sma, e, inc, ape, ran, tru, m, ps)
    dQmin   = Inf
    dQmax   = -Inf
    dQ      = dQn(sma,e,inc,ape,ran,tru,m,ps)
    θs      = range(0.0,2*pi; length = ps.steps)
    for i in eachindex(θs)
        dQθ = dQn(sma,e,inc,ape,ran,θs[i],m,ps)
        if dQθ < dQmin; dQmin = dQθ; end
        if dQθ > dQmax; dQmax = dQθ; end
    end
    ηr      = (dQ - dQmax) / (dQmin - dQmax)
end