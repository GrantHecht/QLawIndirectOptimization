
function qLawEOMsSundmanTransformedZOH(u,p,L)
    # Get time and mass
    t       = u[6]
    m       = u[7]

    # Construct MEE state
    mee     = SVector(u[1], u[2], u[3], u[4], u[5], L)

    # Compute acceleration due to thrust
    if p.coasting == true || p.eclipsed == true
        at  = SVector(0.0, 0.0, 0.0)
        am  = 0.0
    else
        α   = p.α
        β   = p.β
        T   = p.T
        at  = SVector(T*cos(β)*sin(α) / m,
                      T*cos(β)*cos(α) / m,
                      T*sin(β) / m)
        am  = T
    end

    # Compute state dynamics
    dmee    = AstroEOMs.meeEomControl(mee,p.meePs,t,at)
    dLinv   = 1.0 / dmee[6]

    # Compute sundman transformed dynamics
    dmees   = SVector(dmee[1]*dLinv, dmee[2]*dLinv, dmee[3]*dLinv,
                dmee[4]*dLinv, dmee[5]*dLinv, dLinv)

    # Compute mass dynamics
    dm      = -am / p.c
    dms     = dm * dLinv

    # Return full state dynamics
    return SVector(dmees[1],dmees[2],dmees[3],dmees[4],dmees[5],dmees[6],dms)
end

# Our dynamics do not need all the info contained within the qLawParams struct
# Therefore, this method employs a simpler parameter tuple 
function qLawEOMsSundmanTransformedZOH(u, p::Tuple, L, mee_params)
    # Get the parameters
    α = p[1]
    β = p[2]
    T = p[3]
    c = p[4]
    coasting_flag = p[5]

    # Get time and mass
    t       = u[6]
    m       = u[7]

    # Construct MEE state
    mee     = SVector(u[1], u[2], u[3], u[4], u[5], L)

    # Compute acceleration due to thrust
    if coasting_flag
        at  = SA[0.0, 0.0, 0.0]
        am  = 0.0
    else
        at  = SA[
            T*cos(β)*sin(α) / m,
            T*cos(β)*cos(α) / m,
            T*sin(β) / m,
        ]
        am  = T
    end

    # Compute state dynamics
    dmee    = AstroEOMs.meeEomControl(mee, mee_params, t, at)
    dLinv   = 1.0 / dmee[6]

    # Compute sundman transformed dynamics
    dmees   = SA[
        dmee[1]*dLinv, dmee[2]*dLinv, dmee[3]*dLinv,
        dmee[4]*dLinv, dmee[5]*dLinv, dLinv,
    ]

    # Compute mass dynamics
    dm      = -am / c
    dms     = dm * dLinv

    # Return full state dynamics
    return SA[dmees[1],dmees[2],dmees[3],dmees[4],dmees[5],dmees[6],dms]
end