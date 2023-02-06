
function quickestDescentSolve(dQdxA,atMax)
    # Create optimization variables
    u   = Convex.Variable(3)

    # Construct constraint
    con = Convex.norm2(u) <= atMax

    # Construct optimization problem
    p   = Convex.minimize(Convex.dot(dQdxA, u), con)

    # Solve the optimization problem
    Convex.solve!(p, SCS.Optimizer; silent_solver = true)

    if p.status == MOI.OPTIMAL || p.status == MOI.ALMOST_OPTIMAL
        usol = SVector(Convex.evaluate(u)...)
    else
        usol = SVector(0.0, 0.0, 0.0)
    end

    return usol
end