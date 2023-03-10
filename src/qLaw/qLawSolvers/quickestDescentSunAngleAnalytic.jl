

function quickestDescentSunAngleAnalyticSolve(u0,dQdxA,atMax,mee,ps)
    # Get vector a
    a   = dQdxA

    if ps.thrustSunAngleConstraint == false
        na   = norm(a)
        usol = SVector(-atMax*a[1] / na,
                       -atMax*a[2] / na,
                       -atMax*a[3] / na)
    else
        # Determine best direction for cone constraint (sun facing or opposing)
        γ   = tan(pi/2.0 - ps.thrustSunAngle)
        uz  = ps.toSunVec / norm(ps.toSunVec)
        if dot(uz, a) > 0.0
            uz .*= -1.0
        end

        # Construct rotation from LVLH to Inertial
        cart    = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Cartesian, ps.μ)
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

        # Rotate a into cone frame
        ac      = Rls*a

        # Construct matricies for cone constraints
        C       = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
        d       = SVector{3}(0.0, 0.0, 1.0 / γ)

        if norm(C*ac) + dot(d,ac) <= 0.0
            na      = norm(ac)
            usolc   = SVector(-atMax*ac[1]/na,
                              -atMax*ac[2]/na,
                              -atMax*ac[3]/na)
            usol    = transpose(Rls)*usolc

        else # Project onto surface of cone
            den     = sqrt((1 + γ^2)*ac[1]^2 + (1 + γ^2)*ac[2]^2)
            usolc   = SVector(-atMax*ac[1]/den,
                              -atMax*ac[2]/den,
                               atMax*γ*sqrt(ac[1]^2 + ac[2]^2)/den)
            usol    = transpose(Rls)*usolc
        end
    end
    return usol
end