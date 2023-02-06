
# Define cost function
steepestDescentCost(dQdx,A,B,u) = transpose(dQdx)*(B + A*u) / norm(B + A*u)

# Define functions for NLopt
function nlOptCost(x,grad,dQdx,A,B)
    if length(grad) > 0
        zgrad = Zygote.gradient(u -> steepestDescentCost(dQdx,A,B,u), x)[1]
        grad[1] = zgrad[1]
        grad[2] = zgrad[2]
        grad[3] = zgrad[3]
    end
    return steepestDescentCost(dQdx, A, B, x)
end
function nlOptConstraint(x,grad,atMax)
    nx = norm(x)
    if length(grad) > 0
        grad[1] = x[1] / nx 
        grad[2] = x[2] / nx
        grad[3] = x[3] / nx
    end
    return nx - atMax
end

# Optimization function
function steepestDescentControl(u0,dQdx,A,B,atMax)
    # Define Optimizer
    opt     = Opt(:LD_SLSQP, 3)
    
    # Set bounds
    opt.lower_bounds = [-atMax, -atMax, -atMax]
    opt.upper_bounds = [ atMax,  atMax,  atMax]
    opt.xtol_rel     = 1e-6

    # Test functions
    nlOptCost(u0,zeros(3),dQdx,A,B)
    nlOptConstraint(u0,zeros(3),atMax)

    # Set objective function
    opt.min_objective = (x,g) -> nlOptCost(x,g,dQdx,A,B)

    # Set inequality constraint
    inequality_constraint!(opt, (x,g) -> nlOptConstraint(x,g,atMax), 1e-8)

    # Optimize
    (minf, minx, ret) = optimize(opt, u0)

    return SVector(minx[1],minx[2],minx[3])
end

# Optimization function
function steepestDescentControl_SCP(u0,dQdx,A,B,atMax; silent = false)
    # Define some trust region parameters
    ρ           = 1e-5*atMax
    α           = 0.5
    β           = 0.95
    tol         = 1e-8
    maxIters    = 100
    maxLSIters  = 100

    # Define variables
    u           = Convex.Variable(3)

    # Define cost function
    J(u)        = steepestDescentCost(dQdx,A,B,u)

    # Begin the optimization loop
    k           = 0
    done        = false
    uk          = zeros(3)
    ukp1a       = zeros(3)
    ukp1        = zeros(3)
    uk         .= u0
    step        = ρ

    if !silent
        print("Iteration #\tstep\n")
    end
    while !done
        # Zeroth order term
        Jk      = J(uk)

        # First order taylor series term
        dJk     = SVector(Zygote.gradient(J, uk)...)

        # Second order taylor series term
        d2Jk    = Zygote.hessian(J, uk)

        # Extract positive definite part with spectral decomposition
        F       = eigen(d2Jk)
        λs      = F.values
        U       = F.vectors
        for i in eachindex(λs)
            if λs[i] < 0.0
                λs[i] = 0.0
            end
        end
        P       = U*diagm(λs)*transpose(U)

        # Force symmetry
        P[2,1]  = P[1,2]
        P[3,1]  = P[1,3]
        P[3,2]  = P[2,3]

        # Check for convergence
        if k != 0
            γ       = 5.0*step # Not sure if this is the appropriate choice 

            p       = Convex.minimize(
                        transpose(dJk)*(u - uk) + 0.5*Convex.quadform(u - uk, P) +
                        0.5*γ*Convex.norm2(u - uk), 
                        Convex.norm2(u) <= atMax)
            Convex.solve!(p, SCS.Optimizer; silent_solver = true)

            ukp1   .= Convex.evaluate(u)
            G       = SVector(γ*(uk[1] - ukp1[1]), γ*(uk[2] - ukp1[2]), γ*(uk[3] - ukp1[3]))
            if norm(G) <= tol
                done = true
                break
            end
        end

        # Construct constraints
        cons    = Convex.Constraint[
                Convex.norm2(u) <= atMax,
                Convex.norm2(u - uk) <= ρ
                ]

        # Setup problem
        p       = Convex.minimize(
                    transpose(dJk)*(u - uk) + 0.5*Convex.quadform(u - uk, P), 
                    cons)

        # Solve problem
        Convex.solve!(p, SCS.Optimizer; silent_solver = true)

        # Perform backtracking line search to fine next iterate
        ukp1a  .= Convex.evaluate(u)
        Δ       = SVector(ukp1a[1] - uk[1], ukp1a[2] - uk[2], ukp1a[3] - uk[3])
        Jdiff   = -Convex.evaluate(transpose(dJk)*(u - uk) + 0.5*Convex.quadform(u - uk, P))
        t       = 1.0
        l       = 0
        stepAccepted = false
        while !stepAccepted && l <= maxLSIters
            l    += 1
            ukp1t = SVector(uk[1] + t*Δ[1], uk[2] + t*Δ[2], uk[3] + t*Δ[3]) 
            Jkp1a = J(ukp1t)
            if Jkp1a <= Jk - α*t*Jdiff
                stepAccepted = true
                ukp1        .= ukp1t
                step         = sqrt((ukp1[1] - uk[1])^2 + (ukp1[2] - uk[2])^2 + (ukp1[3] - uk[3])^2)
            else
                t   *= β
            end
        end
        uk  .= ukp1
        k   += 1

        if k >= maxIters
            done = true
        end

        # Print status
        if !silent
            print(string(k) * "\t\t" * 
              string(round(step; sigdigits = 5)) * "\n")
        end
    end

    return SVector(uk[1], uk[2], uk[2])
end