
# For mayer cost
function adjointControlTransformation(mee,λv,dS,α,dα,β,dβ,λm,pTuple)
    # Compute cartesian state
    cart = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Cartesian, pTuple[2].μ)

    # Compute thrust direction and first time derivative in LVLH frame
    up  = SVector(cos(α)*cos(β),
                  sin(α)*cos(β),
                  sin(β)) 
    dup = SVector(-dα*sin(α)*cos(β) - dβ*cos(α)*sin(β),
                   dα*cos(α)*cos(β) - dβ*sin(α)*sin(β),
                   dβ*cos(β))

    # Compute rotation matrix and its time derivative
    R, dR   = rotationMatrix(cart, pTuple[2]) 

    # Compute thrust direction and unit vector in inertial frame
    u   = R*up
    du  = R*dup + dR*up

    # Compute velocity costate
    λvv = λv * u

    # Compute requirements for position costate
    c   = pTuple[1]*pTuple[2].TU / (pTuple[2].LU * 1000.0)
    T   = pTuple[1].tMax*pTuple[2].TU^2 / (pTuple[2].MU * pTuple[2].LU * 1000.0) 
    m   = mee[7]
    S   = (c / m)*λv + λm 
    u   = 0.5*(1.0 + tanh(S / pTuple[3].ϵ))
    dm  = -(T / m)*u
    dλm = -(T / m^2)*λv*u
    dλv = (dm / c)*(S - λm) + (m / c)*(dS - dλm)

    # Compute position costate
    λrv = dλv*u + λv*du
    λcart = vcat(λrv, λvv)

    # Convert cartesian co-states to MEE
    dcartdmee   = AstroUtils.convertStatePartials(mee, AstroUtils.MEE, AstroUtils.Cartesian, pTuple[2].μ)
    λmee        = transpose(dcartdmee)*λcart 

    # Return MEE costates
    return vcat(λmee, λm)
end