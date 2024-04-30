
function get_lvlh_to_inertial_mee(mee, ps)
    cart = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Cartesian, ps.μ)
    r    = cart[SA[1,2,3]]
    v    = cart[SA[4,5,6]]

    nr   = norm(r)
    ur   = SA[r[1] /nr, r[2] / nr, r[3] / nr]

    h    = cross(r,v)
    nh   = norm(h)
    uh   = SA[h[1] / nh, h[2] / nh, h[3] / nh]

    t    = cross(uh, ur)
    nt   = norm(t)
    ut   = SA[t[1] / nt, t[2] / nt, t[3] / nt]

    return SA[ur[1] ut[1] uh[1]; ur[2] ut[2] uh[2]; ur[3] ut[3] uh[3]]
end

function get_lvlh_to_inertial_kep(kep, ps)
    cart = AstroUtils.convertState(kep, AstroUtils.Keplerian, AstroUtils.Cartesian, ps.μ)
    r    = cart[SA[1,2,3]]
    v    = cart[SA[4,5,6]]

    nr   = norm(r)
    ur   = SA[r[1] /nr, r[2] / nr, r[3] / nr]

    h    = cross(r,v)
    nh   = norm(h)
    uh   = SA[h[1] / nh, h[2] / nh, h[3] / nh]

    t    = cross(uh, ur)
    nt   = norm(t)
    ut   = SA[t[1] / nt, t[2] / nt, t[3] / nt]

    return SA[ur[1] ut[1] uh[1]; ur[2] ut[2] uh[2]; ur[3] ut[3] uh[3]]
end

function qLaw_control(sma, e, inc, ran, ape, tru, m, constrained, ps::qLawParams)     
    # Compute max thrust acceleration
    atMax   = ps.tMax / m

    # Compute layapunov fuction partials
    dQdx    = Qpartials_keplerian(sma, e, inc, ran, ape, m, ps)

    # Compute requirements
    A, B    = gaussVarKeplerian(sma, e, inc, ran, ape, tru, ps) 

    # Compute product of dQdx^T * A
    a       = SA[
        dQdx[1]*A[1,1] + dQdx[2]*A[2,1] + dQdx[3]*A[3,1] + dQdx[4]*A[4,1] + dQdx[5]*A[5,1],
        dQdx[1]*A[1,2] + dQdx[2]*A[2,2] + dQdx[3]*A[3,2] + dQdx[4]*A[4,2] + dQdx[5]*A[5,2],
        dQdx[1]*A[1,3] + dQdx[2]*A[2,3] + dQdx[3]*A[3,3] + dQdx[4]*A[4,3] + dQdx[5]*A[5,3],
    ]

    # Rotate a into inertial frame
    # Get rotation from lvlh to inertial frame
    RlI = get_lvlh_to_inertial_kep(SA[sma,e,inc,ran,ape,tru], ps)

    # Rotate a into inertial frame
    aI = RlI*a

    # If unconstrained, directly compute
    if !constrained
        neg_a = -a

        # Thrust direction
        a_mag = norm(neg_a)
        dir = SA[neg_a[1] / a_mag, neg_a[2] / a_mag, neg_a[3] / a_mag]
    else
        # Determine best direction for cone constraint (sun facing or opposing)
        γ   = tan(pi/2.0 - ps.thrustSunAngle)
        nts = norm(ps.toSunVec)
        uz  = SA[ps.toSunVec[1] / nts, ps.toSunVec[2] / nts, ps.toSunVec[3] / nts]
        if dot(uz, aI) > 0.0
            uz = -uz
        end

        # Construct rotation matrix to cone-frame
        x  = cross(uz, SA[0,0,1])
        y  = cross(uz, x)
        nx = norm(x)
        ny = norm(y)
        ux = SA[x[1]/nx, x[2]/nx, x[3]/nx]
        uy = SA[y[1]/ny, y[2]/ny, y[3]/ny]
        RIs = SA[ux[1] ux[2] ux[3]; uy[1] uy[2] uy[3]; uz[1] uz[2] uz[3]]

        # Construct full rotation matrix
        Rls = RIs*RlI

        # Rotate a into cone frame
        ac = Rls*a

        # Construct cone constraint matricies
        R = SA[1 0 0; 0 1 0; 0 0 0]
        d = SA[0.0, 0.0, 1.0 / γ]

        if norm(R*ac) + dot(d,ac) <= 0.0
            na      = norm(ac)
            dirc    = SA[-ac[1]/na, -ac[2]/na, -ac[3]/na]
        else # Project onto surface of cone
            den     = sqrt((1 + γ^2)*ac[1]^2 + (1 + γ^2)*ac[2]^2)
            dirc    = SA[-ac[1]/den, -ac[2]/den, γ*sqrt(ac[1]^2 + ac[2]^2)/den]
        end
        dir     = transpose(Rls)*dirc
    end

    # Compute rate of change of Q``
    at      = dir*atMax
    dQdt = transpose(dQdx)*(B + A*at)

    return dir, dQdt
end

function effectivity_check(dQdt, sma, e, inc, ran, ape, m, constrained, ps::qLawParams)
    if ps.ηa > 0.0 || ps.ηr > 0.0
        # Compute dQmin and dQmax
        θs      = range(0.0, 2*pi; length = ps.steps)
        dQmin   = Inf
        dQmax   = -Inf
        for i in eachindex(θs)
            θ = θs[i]
            _, dQdt_θ = qLaw_control(sma, e, inc, ran, ape, θ, m, constrained, ps)
            if dQdt_θ < dQmin; dQmin = dQdt_θ; end
            if dQdt_θ > dQmax; dQmax = dQdt_θ; end
        end

        # Compute effectivity terms
        ηa      = dQdt / dQmin
        ηr      = (dQdt - dQmax) / (dQmin - dQmax)
        ηa_val  = ηa - ps.ηa 
        ηr_val  = ηr - ps.ηr

        if ηa_val < 0.0 || ηr_val < 0.0
            return true
        else
            return false
        end
    else
        return false
    end
end

function qLaw_control_with_coast_checks(mee, m, constrained::Bool, ps::qLawParams)
    # Convert to Keplerian elements
    kep     = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Keplerian, ps.μ)
    sma     = kep[1]
    e       = kep[2]
    inc     = kep[3]
    ran     = kep[4]
    ape     = kep[5]
    tru     = kep[6]

    # Compute thrust direction
    dir, dQdt = qLaw_control(sma, e, inc, ran, ape, tru, m, constrained, ps)

    # Check if our control reduces Q
    positive_qdot_coast = dQdt > 0.0

    # Effectivity check
    effectivity_coast = if positive_qdot_coast
        false
    else
        effectivity_check(dQdt, sma, e, inc, ran, ape, m, constrained, ps)
    end

    # Compute angles
    coast_flag = effectivity_coast || positive_qdot_coast
    if !coast_flag
        α = atan(dir[1],dir[2])
        β = atan(dir[3] / sqrt(dir[1]*dir[1] + dir[2]*dir[2]))
    else
        α = 0.0
        β = 0.0
    end

    return (α,β,ps.tMax,effectivity_coast,positive_qdot_coast)
end
