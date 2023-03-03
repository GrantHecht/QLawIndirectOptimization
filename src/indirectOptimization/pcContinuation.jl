function computeTangentVector(J)
    # Compute nullspace of J
    nullVec = nullspace(J)

    # Check determinate condition
    if det(vcat(J, transpose(nullVec))) > 0.0
        return nullVec
    else 
        return -nullVec
    end
end


function pcContinuation(λ0, Δs, mee0, interp, pTuple)
    # Preallocate for forward diff
    result1 = DiffResults.JacobianResult(zeros(6), rand(7))
    mem1    = zeros(6)
    config1 = ForwardDiff.JacobianConfig((F,x) -> minFuelMayerNonlinearFunctionWithInterp!(F,x,mee0,interp,pTuple),
                zeros(6), rand(7))
    result2 = DiffResults.JacobianResults(zeros(7), rand(7))
    min2    = zeros(7)
    config2 = ForwardDiff.JacobianConfig((F,x) -> minFuelMayerPseudoarclengthNonlinearFunction!(F,x,zeros(7),zeros(7),Δs,mee0,interp,pTuple),
                zeros(7), rand(7))

    # Allocate iteration variables
    u      = zeros(7)
    v      = zeros(7)
    F      = zeros(6)
    J      = zeros(6, 7)

    # Begin loop
    u[1:6] .= λ0
    u[7]    = 0.0
    h       = h0
    i       = 0
    while i < 1
        # Increment iteration counter
        i += 1

        # Get Jacobian of nonlinear function at current step 
        result1 = ForwardDiff.jacobian!(result1, (F,x) -> minFuelMayerNonlinearFunctionWithInterp!(F,x,mee0,interp,pTuple),
                                        mem1, u, config1, Val{false}())
        F      .= DiffResults.value(result1)
        J      .= DiffResults.jacobian(result1)

        # Compute tangent vector of J
        dc      = computeTangentVector(J)

        # Compute initial guess
        v      .= u .+ Δs*dc

        # Solve problem with NLsolve
        sol     = nlsolve(only_fj!(F,J,x) -> minFuelMayerPseudoarclengthNonlinearFunctionWithJacobian!(F,J,x,u,dc,Δs,mee0,meef_interp,pTuple,result2,mem2,config2),
                    v; show_trace = true)
    end
end