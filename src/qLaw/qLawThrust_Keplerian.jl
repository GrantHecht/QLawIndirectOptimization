
function qLawThrust_Keplerian(mee, m, ps::qLawParams)
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
    p       = sma*(1.0 - e*e)
    h       = sqrt(ps.μ*p)
    r       = p / (1.0 + e*cos(tru))
    hinv    = 1.0 / h

    # Compute variational equations (without ta)
    A       = SMatrix{6,3}(2*sma*sma*hinv*e*sin(tru),       # dadfr
                           hinv*p*sin(tru),                 # dedfr
                           0.0,                             # didfr
                           0.0,                             # dΩdfr
                           -p*hinv*cos(tru) / e,            # dΩdfr
                           p*hinv*cos(tru) / e,             # dθdfr
                           2*sma*sma*hinv*p / r,            # dadfθ
                           hinv*((p + r)*cos(tru) + r*e),   # dedfθ
                           0.0,                             # didfθ 
                           0.0,                             # dΩdfθ
                           hinv*(p + r)*sin(tru) / e,       # dωdfθ
                           -(p + r)*hinv*sin(tru) / e,      # dθdfθ
                           0.0,                             # dadfh
                           0.0,                             # dedfh
                           r*hinv*cos(tru + ape),           # didfh
                           r*hinv*sin(tru + ape) / sin(inc),# dΩdfh
                           -(r*hinv*sin(tru + ape)*cos(inc) / sin(inc)),
                           0.0)

    B       = SVector(0.0, 0.0, 0.0, 0.0, 0.0, h / (r*r))

    # Compute product of dQdx^T * A
    dQdxA   = SVector(dQdx[1]*A[1,1] + dQdx[2]*A[2,1] + dQdx[3]*A[3,1] + dQdx[4]*A[4,1] + dQdx[5]*A[5,1],
                      dQdx[1]*A[1,2] + dQdx[2]*A[2,2] + dQdx[3]*A[3,2] + dQdx[4]*A[4,2] + dQdx[5]*A[5,2],
                      dQdx[1]*A[1,3] + dQdx[2]*A[2,3] + dQdx[3]*A[3,3] + dQdx[4]*A[4,3] + dQdx[5]*A[5,3])

    # Solve the steapest decent optimization problem first
    #   Minimize:   J = dQdxA * u 
    #   Subject to: ||u|| <= atMax

    # Variables
    u   = Convex.Variable(3)

    # Construct problem with cost
    p   = Convex.minimize(Convex.dot(dQdxA, u))

    # Add constraint
    p.constraints += Convex.norm2(u) <= atMax

    # Solve
    Convex.solve!(p, SCS.Optimizer; silent_solver = true)

    # If solve is successfull, move to sequential quadratic optimization loop
    if p.status == MOI.OPTIMAL || p.status == MOI.ALMOST_OPTIMAL
        uk  = Convex.evaluate(u)

        # Define nonlinear cost function
        dQdxf       = SVector(dQdx[1], dQdx[2], dQdx[3], dQdx[4], dQdx[5], 0.0)
        J(u)        = transpose(dQdxf)*(B + A*u) / norm(B + A*u)

        # Define some trust region parameters
        ρ           = 0.00001 * atMax
        tol         = 1e-6
        maxiters    = 50

        # Begin optimization loop
        k    = 0
        ukp1 = zeros(3)
        done = false

        print("Iteration #    diff\n")
        while !done
            # Zeroth order taylor series term
            Jk     = J(uk)
            
            # First order taylor series term
            dJk   = SVector(Zygote.gradient(J, uk)...)

            # Second order taylor series term
            d2Jk   = Zygote.hessian(J, uk)

            # Compute eigen decomposition of second order term
            F      = eigen(d2Jk)
            λs     = F.values
            U      = F.vectors
            for i in eachindex(λs)
                λs[i] = λs[i] < 0.0 ? 0.0 : λs[i]
            end
            display(λs)
            P      = U * diagm(λs) * transpose(U)
            display(P)

            # Force symmetry
            P[2,1] = P[1,2]
            P[3,1] = P[1,3]
            P[3,2] = P[2,3]

            # Construct constraints
            cons   = Convex.Constraint[
                Convex.norm2(u) <= atMax,
                Convex.norm2(u - uk) <= ρ
            ]

            # Setup problem
            p      = Convex.minimize(transpose(dJk)*(u - uk) + 0.5*Convex.quadform(u - uk, P), cons)
            #p      = Convex.minimize(transpose(dJk)*(u - uk), cons)

            # Solve the problem
            Convex.solve!(p, SCS.Optimizer; silent_solver = true, warmstart = k == 0 ? false : true)

            # Get solution
            ukp1 .= Convex.evaluate(u)

            # Compute differene
            diff = norm(ukp1 - uk) / norm(uk)

            # Update uk
            uk .= ukp1

            # Increment iteration counver
            k += 1

            # Print status
            print(string(k) * "          " * string(diff) * "\n")

            # Checking for stopping conditions
            if diff < tol
                done = true
            elseif k >= maxiters
                done = true
            end
        end

        uopt    = SVector(m*uk[1], m*uk[2], m*uk[3])
    else
        uopt    = SVector(0.0, 0.0, 0.0)
    end

    #return (p.status, uopt)
    return uopt
end