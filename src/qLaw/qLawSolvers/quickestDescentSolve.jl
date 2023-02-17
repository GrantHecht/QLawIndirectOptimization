wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
function quickestDescentSolve(u0,dQdxA,atMax,ps)
    # Create optimization variables
    u   = Convex.Variable(3)
    u.value = u0

    # Construct constraint
    if ps.thrustSunAngleConstraint == false
        cons = [Convex.norm2(u) <= atMax]

        # Construct optimization problem
        p   = Convex.minimize(Convex.dot(dQdxA, u), cons)

        #@suppress begin
        Convex.solve!(p, SCS.Optimizer; silent_solver = false, warmstart = false)
        #end
        wait_for_key("First cone solve finished")

        if p.status == MOI.OPTIMAL || p.status == MOI.ALMOST_OPTIMAL
            usol = SVector(Convex.evaluate(u)...)
        else
            usol = SVector(0.0, 0.0, 0.0)
        end
    else
        # Determine best direction for cone constraint (sun facing or opposing)
        a   = tan(pi/2.0 - ps.thrustSunAngle)
        uz  = ps.toSunVec / norm(ps.toSunVec)
        if dot(uz, dQdxA) > 0.0
            uz .*= -1.0
        end

        # Construct rotation matrix to sun-frames (improve for efficiency later)
        ux = cross(uz, [0, 0, 1])
        uy = cross(uz,ux)
        A  = SMatrix{3,3}(ux[1], uy[1], 0, 
                          ux[2], uy[2], 0, 
                          ux[3], uy[3], 0)
        d  = SVector(uz[1] / a,uz[2] / a,uz[3] / a)

        # Constraints for sun facing cone
        cons = [Convex.norm2(u) <= atMax,
                Convex.norm2(A*u) <= Convex.dot(d,u),
                Convex.dot(d,u) >= 0.0]

        # Problem
        p = Convex.minimize(Convex.dot(dQdxA, u), cons)

        @suppress begin
        Convex.solve!(p, SCS.Optimizer; silent_solver = false, warmstart = false)
        end

        # If problem status is dual infeasible, try with alternative cone constraint
        if p.status == MOI.DUAL_INFEASIBLE
            # Construct rotation matrix to sun-frames (improve for efficiency later)
            ux .*= -1.0
            ux = cross(uz, [0, 0, 1])
            uy = cross(uz,ux)
            A  = SMatrix{3,3}(ux[1], uy[1], 0, 
                            ux[2], uy[2], 0, 
                            ux[3], uy[3], 0)
            d  = SVector(uz[1] / a,uz[2] / a,uz[3] / a)

            # Constraints for sun facing cone
            cons = [Convex.norm2(u) <= atMax,
                    Convex.norm2(A*u) <= Convex.dot(d,u),
                    Convex.dot(d,u) >= 0.0]

            # Problem
            p = Convex.minimize(Convex.dot(dQdxA, u), cons)

            @suppress begin
            Convex.solve!(p, SCS.Optimizer; silent_solver = false, warmstart = false)
            end
        end

        if p.status == MOI.OPTIMAL || p.status == MOI.ALMOST_OPTIMAL
            usol = SVector(Convex.evaluate(u)...)
        else
            usol = SVector(0.0, 0.0, 0.0)
        end
    end

    return usol
end