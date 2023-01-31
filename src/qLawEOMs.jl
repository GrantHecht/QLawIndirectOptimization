
function qLawEOMs(u, p, L)
    # Get states
    p   = u[1]
    e   = u[2]
    inc = u[3]
    ran = u[4]
    ape = u[5]
    m   = u[6]
    t   = u[7]

    # Other orbital elements
    tru = L - ran - ape
    sma = p / (1.0 - e*e)
    h   = sqrt(ps.μ*p)
    r   = p / (1.0 + e*cos(tru))
    g   = e*sin(ape + ran)
    w   = 1.0 + e*cos(ape + ran)*cos(L) * g*sin(L)
    k   = tan(inc / 2.0)*sin(ran)

    # Thrust and acceleration
    current_thrust = ps.tMax
    f   = current_thrust / m
    
    # Compute Lyapunov Control
    fr, fθ, fh = qLawThrust(sma, e, inc, ape, ran, tru, m, ps)
    nf         = sqrt(fr*fr + fθ*fθ + fh*fh)
    fr         /= nf
    fθ         /= nf
    fh         /= nf

    # Effectivity calculation
    coast = false
    if ps.effectivity

    end
end