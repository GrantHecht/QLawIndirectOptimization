
function qLawEffectivity_Keplerian(mee, m, ps::qLawParams; method = :SD)
    # Convert to Keplerian elements
    kep,f   = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Keplerian, ps.μ)
    sma     = kep[1]
    e       = kep[2]
    inc     = kep[3]
    ran     = kep[4]
    ape     = kep[5]
    tru     = kep[6]

    # Compute max thrust acceleration
    atMax   = ps.tMax / m

    # Compute layapunov fuction partials
    dQdx    = Qpartials_keplerian(sma, e, inc, ran, ape, m, ps)

    # Compute dQ terms
    dQmin   = Inf
    dQmax   = -Inf
    θs      = range(0.0, 2*pi; length = ps.steps)
    for i in eachindex(θs)
        # Compute gauss variational equations
        (A,B)   = gaussVarKeplerian(sma, e, inc, ran, ape, θs[i], ps)

        # Compute product of dQdx^T * A
        dQdxA   = SVector(dQdx[1]*A[1,1] + dQdx[2]*A[2,1] + dQdx[3]*A[3,1] + dQdx[4]*A[4,1] + dQdx[5]*A[5,1],
                        dQdx[1]*A[1,2] + dQdx[2]*A[2,2] + dQdx[3]*A[3,2] + dQdx[4]*A[4,2] + dQdx[5]*A[5,2],
                        dQdx[1]*A[1,3] + dQdx[2]*A[2,3] + dQdx[3]*A[3,3] + dQdx[4]*A[4,3] + dQdx[5]*A[5,3])

        # If using quickest decent with no constraint, directly compute
        if method == :QDUC
            # Unconstrianed thrust 
            atQDUC = -dQdxA

            # Thrust direction
            aMag = norm(atQDUC)
            dir  = SVector(atQDUC[1] / aMag, atQDUC[2] / aMag, atQDUC[3] / aMag)

            # Thrust magnitude
            tMag = ps.tMax
        else
            # Solve the quickest decent optimization problem first
            #atQD = quickestDescentSolve(dQdxA, atMax)
            atQD = -dQdxA

            # Solve with steepest descent using quickest decent solution as guess
            if method == :SD
                # Ensure initial guess is in bounds
                mag = norm(atQD)
                dir = SVector(atQD[1] / mag, atQD[2] / mag, atQD[3] / mag)
                u0  = SVector(atMax*dir[1], 
                            atMax*dir[2], 
                            atMax*dir[3])

                atSD = steepestDescentControl(u0,dQdx,A,B,atMax)

                # Check that control is reducing Q
                dQdx_dot_fn = transpose(dQdx)*(B + A*atSD)/norm(B + A*atSD)
                if dQdx_dot_fn > 0.0
                    atSD .= 0.0
                end

                # COmpute thrust direction
                aMag = norm(atSD)
                dir  = SVector(atSD[1]/aMag, atSD[2]/aMag, atSD[3]/aMag)

                # Compute thrust magnitude
                tMag = m*aMag
                if tMag > ps.tMax
                    tMag = ps.tMax
                end
            else
                # Check that control is reducing Q
                dQdt = transpose(dQdx)*(B + A*atQD)
                if dQdt > 0.0
                    atQD .= 0.0
                end

                # Compute thrust direction
                aMag = norm(atQD)
                dir  = SVector(atQD[1]/aMag, atQD[2]/aMag, atQD[3]/aMag)

                # Compute thrust magnitude
                tMag = m*aMag
                if tMag > ps.tMax
                    tMag = ps.tMax
                end
            end
        end

        at      = SVector(tMag*dir[1]/m, tMag*dir[2]/m, tMag*dir[3]/m)
        f       = B .+ A*at
        if method != :SD
            dQt     = transpose(dQdx)*f
        else
            dQt     = transpose(dQdx)*f/norm(f)
        end
        if dQt < dQmin; dQmin = dQt; end
        if dQt > dQmax; dQmax = dQt; end
    end

    return (dQmin, dQmax)
end