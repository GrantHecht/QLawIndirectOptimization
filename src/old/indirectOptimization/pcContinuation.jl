function computeTangentVector(J, posDetDirection)
    # Compute nullspace of J
    nulls   = nullspace(J)

    # Grab nullsapce vector with smallest product norm if required
    nullVec = zeros(size(J, 2))
    if size(nulls, 2) > 1
        midx    = 0
        mnorm   = Inf
        prod    = zeros(size(J, 1))
        @inbounds for i in axes(nulls, 2)
            mul!(prod, J, view(nulls, :, i))
            if norm(prod) < mnorm
                midx = i
            end
            nullVec .= nulls[:,midx]
        end
    else
        nullVec .= nulls
    end


    # Check determinate condition
    if det(vcat(J, transpose(nullVec))) > 0.0
        if posDetDirection
            return nullVec
        else
            return -nullVec
        end
    else 
        if posDetDirection
            return -nullVec
        else
            return nullVec
        end
    end
end

function pcContinuation(λ0, Δs0, mee0, interp, pTuple; critPointTol = 1e-16)
    # Preallocate for forward diff
    result1 = DiffResults.JacobianResult(zeros(6), rand(7))
    mem1    = zeros(6)
    config1 = ForwardDiff.JacobianConfig((F,x) -> minFuelMayerNonlinearFunctionWithInterp!(F,x,mee0,interp,pTuple),
                zeros(6), rand(7))
    result2 = DiffResults.JacobianResult(zeros(7), rand(7))
    mem2    = zeros(7)
    config2 = ForwardDiff.JacobianConfig((F,x) -> minFuelMayerPseudoarclengthNonlinearFunction!(F,x,zeros(7),zeros(7),Δs,mee0,interp,pTuple),
                zeros(7), rand(7))

    # Allocate iteration variables
    u      = zeros(7)
    v      = zeros(7)
    F      = zeros(6)
    J      = zeros(6, 7)

    # Allocate variables for debugging
    κs     = zeros(1000)

    # Begin loop
    Δs          = Δs0
    u[1:6]     .= λ0
    u[7]        = 0.0
    i           = 0
    checkDetDir = true
    posDetDir   = true
    while i < 1000
        # Increment iteration counter
        i += 1

        # Get Jacobian of nonlinear function at current step 
        result1 = ForwardDiff.jacobian!(result1, (F,x) -> minFuelMayerNonlinearFunctionWithInterp!(F,x,mee0,interp,pTuple),
                                        mem1, u, config1, Val{false}())
        F      .= DiffResults.value(result1)
        J      .= DiffResults.jacobian(result1)

        # Compute determinate of shooting function jacobian
        #detJF   = det(view(J, 1:6, 1:6))
        rcondJ  = 1.0 / cond(view(J, 1:6, 1:6), 1)

        # If determinate is less than critPointTol, perform fixed-point homotopy to find new path
        println("rcondJ = " * string(rcondJ))
        println("Δs = " * string(Δs))
        if rcondJ < critPointTol
            println("Beginning fixed point homotopy.")
            newSolutionAccepted = false
            while !newSolutionAccepted
                v      .= fpContinuation(u, 1e-2, mee0, interp, pTuple)
                v[end]  = u[end]

                newSolutionAccepted = true
                checkDetDir         = true
                u                  .= v

                # Solve with new u (May want to add corection at some point, but I think the below code needs to be modified)
                #sol     = nlsolve(only_fj!((F,J,x) -> minFuelMayerPseudoarclengthNonlinearFunctionWithJacobian!(F,J,x,u,zeros(7),0.0,mee0,interp,pTuple,result2,mem2,config2)),
                #            v; show_trace = false)

                #if sol.f_converged == true
                #    println("New solution accepted.")
                #    newSolutionAccepted = true
                #    #u .= v
                #    u .= sol.zero
                #end
            end
        else # Continue along current path
            # Compute tangent vector of J
            dc      = computeTangentVector(J, posDetDir)
            
            # If necessary, make sure we're stepping with increasing κ1
            if checkDetDir 
                checkDetDir = false
                if u[7] + Δs*dc[7] < u[7]
                    posDetDir = !posDetDir
                    dc  = computeTangentVector(J, posDetDir)
                end
            end

            # Compute initial guess
            v      .= u .+ Δs*dc

            # Solve problem with NLsolve
            sol     = nlsolve(only_fj!((F,J,x) -> minFuelMayerPseudoarclengthNonlinearFunctionWithJacobian!(F,J,x,u,dc,Δs,mee0,interp,pTuple,result2,mem2,config2)),
                        v; show_trace = false)

            # If successful, update solution
            if sol.f_converged == true
                Δs *= 1.1
                κs[i] = sol.zero[end]
                println("κ1 = " * string(sol.zero[end]))
                u  .= sol.zero
            else
                println("Reducing Δs!")
                Δs /= 2.0
            end
        end
    end
end

function fpContinuation(λ0, Δs, mee0, interp, pTuple; exitInitialIterateTol = 1e-7, solutionTol = 1e-12)
    # Grab componants of λ0
    y1  = λ0[1:6]
    κ1  = λ0[7]

    # Iniitialize second continuation term
    κ2  = 1.0

    # Preallocate for ForwardDiff
    result  = DiffResults.JacobianResult(zeros(7), rand(7))
    mem     = zeros(7)
    config  = ForwardDiff.JacobianConfig((F,x) -> minFuelMayerFixedPointNonlinearFunction!(F,x,zeros(7),zeros(7),Δs,y1,κ1,mee0,interp,pTuple),
                zeros(7), rand(7))

    # Allocate iteration variables
    u      = zeros(7)
    v      = zeros(7)
    F      = zeros(7)
    J      = zeros(7, 7)

    # Begin loop
    u[1:6]     .= y1
    u[7]        = κ2
    i           = 0
    done        = false
    maxIters    = 10000
    hasExitedInitialIterate = false
    while i < maxIters && !done
        # Increment iteration counter
        i += 1

        # Get Jacobian of nonlinear function at current step 
        result  = ForwardDiff.jacobian!(result, (F,x) -> minFuelMayerFixedPointNonlinearFunction!(F,x,u,u,Δs,y1,κ1,mee0,interp,pTuple),
                                        mem, u, config, Val{false}())
        F      .= DiffResults.value(result)
        J      .= DiffResults.jacobian(result)

        # Compute determinate of shooting function jacobian
        detJF   = det(view(J, 1:6, 1:6))

        # If determinate is less than critPointTol, perform fixed-point homotopy to find new path
        #println("detJ = " * string(detJF))

        # Compute tangent vector of J
        dc      = computeTangentVector(view(J, 1:6, :), true)

        # Compute initial guess
        if hasExitedInitialIterate == false
            # Standard initial step
            v      .= u .+ Δs*dc
        else
            # Adjust initial step so we don't step past new solution
            if (u[7] > 1.0 && u[7] + Δs*dc[7] < 1.0) || (u[7] < 1.0 && u[7] * Δs*dc[7] > 1.0)
                Δs = (1.0 - u[7]) / dc[7]
                println("New Δs = " * string(Δs))
                v .= u .+ Δs*dc
            else 
                v .= u .+ Δs*dc
            end
        end

        # Solve problem with NLsolve
        sol     = nlsolve(only_fj!((F,J,x) -> minFuelMayerFixedPointNonlinearFunctionWithJacobian!(F,J,x,u,dc,Δs,y1,κ1,mee0,interp,pTuple,result,mem,config)),
                    v; show_trace = false)

        # If successful, update solution
        if sol.f_converged == true
            Δs *= 1.1
            #println("κ2 = " * string(sol.zero[end]))
            u  .= sol.zero

            # If beginning, check if we've left initial iterate
            if hasExitedInitialIterate == false
                if abs(u[7] - 1.0) > exitInitialIterateTol
                    println("Have exited initial iterate.")
                    hasExitedInitialIterate = true
                end

            # Otherwise, see if we've returned to 1.0
            elseif abs(u[7] - 1.0) < solutionTol
                done = true
            end
        else
            println("Reducing Δs!")
            Δs /= 2.0
        end
    end
    return u
end