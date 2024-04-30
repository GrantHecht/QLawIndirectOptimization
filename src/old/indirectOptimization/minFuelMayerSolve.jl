# Solve minimum fuel problem
function minFuelMayerSolve(mee0, meef, tspan, pTuple; 
                           λ0g = nothing, ϵs = 1.0, 
                           attempts = 100, show_trace = false, 
                           ftol = 1e-8, iterations = 1000)
    # Configure inputs
    result  = DiffResults.JacobianResult(rand(6))
    mem     = zeros(6)
    config  = ForwardDiff.JacobianConfig((F,x) -> minFuelMayerNonlinearFunction!(F,x,mee0,meef,tspan,pTuple), 
                                         zeros(6), rand(6))

    λ0      = zeros(6)
    if !isnothing(λ0g) # If guess was given, use to solve
        λ0 .= λ0g
        λ0s, retcode = minFuelMayerSolve(λ0, ϵs, mee0, meef, tspan, pTuple, result, mem, config, ftol, iterations, show_trace)

    else # Otherwise, randomly guess until convergence is successfull
        count = 0  
        ϵ     = 1.0
        while count < attempts
            # Update counter
            count += 1

            # Generate random guess
            λ0 .= 0.1*randn(6)

            # Attempt to solve
            λ0s, retcode = minFuelMayerSolve(λ0, ϵ, mee0, meef, tspan, pTuple, result, mem, config, ftol, iterations, show_trace)

            # If successfull, break from loop
            if retcode == :success
                break
            end
        end

        # If solve was successful, continue with homotopy continuation
        if retcode == :success && ϵ != ϵs
            λ0s, retcode = minFuelMayerSolve(λ0s, ϵs, mee0, meef, tspan, pTuple, result, mem, config, ftol, iterations, show_trace)
        end
    end

    # Return 
    return (λ0s, retcode)
end

function minFuelMayerSolve(λ0, ϵ, mee0, meef, tspan, pTuple, result, mem, config, ftol, iterations, show_trace)
    # Configure function
    fun     = only_fj!((F,J,x) -> minFuelMayerNonlinearFunctionWithJacobian!(F, J, x, mee0, meef, tspan,
                                                                             pTuple, result, mem, config))

    # Solve
    # Set homotopy parameter
    pTuple[3].ϵ = ϵ

    # Solve 
    nsol    = nlsolve(fun, λ0, show_trace = show_trace, ftol = ftol, iterations = iterations)
    λ0     .= nsol.zero

    # Check if we converged
    if nsol.f_converged
        retcode = :success
    else
        retcode = :failed 
    end

    # Return 
    return (λ0, retcode)
end

function minFuelMayerSolve(λ0, ϵs::AbstractArray, mee0, meef, tspan, pTuple, result, mem, config, ftol, iterations, show_trace)
    # Configure function
    fun     = only_fj!((F,J,x) -> minFuelMayerNonlinearFunctionWithJacobian!(F, J, x, mee0, meef, tspan,
                                                                             pTuple, result, mem, config))

    # Solve
    retcode = :success
    for ϵ in ϵs 
        # Set homotopy parameter
        pTuple[3].ϵ = ϵ

        # Solve 
        nsol    = nlsolve(fun, λ0, show_trace = show_trace, ftol = ftol, iterations = iterations)
        λ0     .= nlsol.zero

        # Check if we converged
        if nsol.f_converged == false
            retcode = :failed 
            break
        end
    end

    # Return 
    return (λ0, retcode)
end