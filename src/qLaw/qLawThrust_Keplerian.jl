
function qLawThrust_Keplerian(mee, m, ps::qLawParams; method = :SD)
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

    # Compute requirements
    A, B    = gaussVarKeplerian(sma, e, inc, ran, ape, tru, ps) 

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
        # Construct guess 
        atQDUC = -dQdxA 
        mag    = norm(atQDUC)
        guess  = SVector(atMax*atQDUC[1] / mag, atMax*atQDUC[2] / mag, atMax*atQDUC[3] / mag)

        # Solve the quickest decent optimization problem first
        if method != :QDSAA
            atQD = quickestDescentSolve(guess, dQdxA, atMax, mee, ps)
        else
            atQD = quickestDescentSunAngleAnalyticSolve(guess, dQdxA, atMax, mee, ps)
        end

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
            if aMag > 0
                dir  = SVector(atSD[1]/aMag, atSD[2]/aMag, atSD[3]/aMag)
            else
                dir  = SVector(0.0, 0.0, 0.0)
            end

            # Compute thrust magnitude
            tMag = m*aMag
            if tMag > ps.tMax
                tMag = ps.tMax
            end
        else
            # Check that control is reducing Q
            dQdt = transpose(dQdx)*(B + A*atQD)
            if dQdt > 0.0
                atQD = SVector(0.0, 0.0, 0.0)
            end

            # Compute thrust direction
            aMag = norm(atQD)
            if aMag > 0
                dir  = SVector(atQD[1]/aMag, atQD[2]/aMag, atQD[3]/aMag)
            else
                dir = SVector(0.0, 0.0, 0.0)
            end

            # Compute thrust magnitude
            tMag = m*aMag
            if tMag > ps.tMax
                tMag = ps.tMax
            end
        end
    end

    # Compute angles
    flag = true
    if norm(dir) > 0.0
        α = atan(dir[1],dir[2])
        β = atan(dir[3] / sqrt(dir[1]*dir[1] + dir[2]*dir[2]))
    else
        flag = false
        α = 0.0
        β = 0.0
    end

    # Effectivity check
    coast = false
    if ps.ηa > 0.0 || ps.ηr > 0.0 && flag == true
        # Compute dQ, dQmin, and dQmax
        at = SVector(dir[1]*tMag/m, dir[2]*tMag/m, dir[3]*tMag/m)
        if method != :SD
            dQ  = transpose(dQdx)*(B .+ A*at)
        else
            f = B .+ A*at
            dQ  = transpose(dQdx)*f / norm(f)
        end
        dQmin, dQmax = qLawEffectivity_Keplerian(mee, m, ps; method = method)

        # Compute effectivity terms
        ηa      = dQ / dQmin
        ηr      = (dQ - dQmax) / (dQmin - dQmax)
        ηa_val  = ηa - ps.ηa 
        ηr_val  = ηr - ps.ηr

        if ηa_val < 0.0 || ηr_val < 0.0
            coast = true
        end
    end

    if flag == false
        coast = true
    end

    return (α,β,tMag,coast)
end