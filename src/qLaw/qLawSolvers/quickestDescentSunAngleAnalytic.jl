

function quickestDescentSunAngleAnalyticSolve(u0,dQdxA,atMax,mee,ps)
    # First compute quickest decent control with only the thrust acceleration constraint
    dQdxAn  = norm(dQdxA)
    usoln   = SVector(-atMax*dQdxA[1] / dQdxAn, -atMax*dQdxA[2] / dQdxAn, -atMax*dQdxA[3] / dQdxAn)

    if ps.thrustSunAngleConstraint == false
        usol = SVector(usoln...)
    else
        # Determine best direction for cone constraint (sun facing or opposing)
        a   = tan(pi/2.0 - ps.thrustSunAngle)
        uz  = ps.toSunVec / norm(ps.toSunVec)
        if dot(uz, dQdxA) > 0.0
            uz .*= -1.0
        end

        # Construct rotation from LVLH to Inertial
        cart    = AstroUtils.convertState(mee,AstroUtils.MEE,AstroUtils.Cartesian,ps.μ)
        ur      = cart[1:3] / norm(view(cart,1:3))
        uh      = cross(view(cart,1:3),view(cart,4:6))
        nh      = norm(uh)
        uh    ./= nh
        ut      = cross(uh, ur)
        RlI     = SMatrix{3,3}(ur[1], ur[2], ur[3],
                               ut[1], ut[2], ut[3],
                               uh[1], uh[2], uh[3])

        # Construct rotation matrix to sun-frames (improve for efficiency later)
        ux      = cross(uz, [0, 0, 1])
        uy      = cross(uz,ux)
        RIs     = SMatrix{3,3}(ux[1], uy[1], uz[1], 
                               ux[2], uy[2], uz[2], 
                               ux[3], uy[3], uz[3])

        # Construct full rotation matrix
        Rls     = RIs*RlI

        # Construct matricies for cone constraints
        A       = SMatrix{2,3}(Rls[1,1], Rls[2,1], Rls[1,2], Rls[2,2], Rls[1,3], Rls[2,3])
        d       = SVector{3}(Rls[3,1] / a, Rls[3,2] / a, Rls[3,3] / a)

        if norm(A*usoln) - transpose(d)*usoln <= 0.0 # Already in cone
            usol = SVector(usoln...)

        else # Project onto surface of cone
            # Rotating control into sun frame
            usolnsun    = Rls*usoln

            # Compute required z componant such that vector is in cone
            uzsun       = a*sqrt(usolnsun[1]^2 + usolnsun[2]^2)

            # Construct new unit vector on surface of cone 
            onconevecn  = sqrt(usolnsun[1]^2 + usolnsun[2]^2 + uzsun^2)
            onconevecu  = SVector(usolnsun[1] / onconevecn, 
                                  usolnsun[2] / onconevecn, 
                                  uzsun / onconevecn)

            # Construct control with correct length
            usolsun     = SVector(atMax*onconevecu[1],
                                  atMax*onconevecu[2],
                                  atMax*onconevecu[3])

            # Rotate back to LVLH
            usol        = transpose(Rls)*usolsun

            # Test
            #println(norm(A*usol) - transpose(d)*usol)
            #wait_for_key("hi")
        end
    end
    return usol
end