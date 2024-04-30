
# Nonlinear function without jacobian
function memfNonlinearFunction!(F, x, mee0, meef, tspan, pTuple)
    # Construct full initial state
    y0  = SVector(mee0[1], mee0[2], mee0[3], mee0[4], mee0[5], mee0[6], pTuple[1].initMass / pTuple[2].MU,
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7])

    # Initialize switching state
    AstroEOMs.meeMEMFInitializeSwitchingState!(pTuple[3], y0, tspan[1], pTuple[1], pTuple[2])

    # Construct callback
    cb   = ContinuousCallback(AstroEOMs.meeMEMFSwitchingCallbackCondition, AstroEOMs.meeMEMFSwitchingCallbackAffect!;
                                save_positions = (false,false),
                                interp_points = 20,
                                abstol = 10eps(),
                                reltol = 0,
                                repeat_nudge = 1//100) 

    # Construct problem and solve
    prob = ODEProblem(ODEFunction{false}(AstroEOMs.meeOptimalControlEOMs), y0, tspan, pTuple; callback = cb)
    sol  = solve(prob, Vern9(), reltol = 1e-14, abstol = 1e-14, save_everystep = false, save_start = false, initialize_save = false, maxiters = 1e7,
            sensealg = ForwardDiffSensitivity(convert_tspan=true))
    xf   = sol[end]

    F[1] = xf[1] - meef[1]
    F[2] = xf[2] - meef[2] 
    F[3] = xf[3] - meef[3]
    F[4] = xf[4] - meef[4]
    F[5] = xf[5] - meef[5]
    #F[6] = xf[6] - meef[6]
    F[6] = xf[13]
    F[7] = xf[14]
end

# Nonlinear function with jacobian
function memfNonlinearFunctionWithJacobian!(F, J, x, mee0, meef, tspan, pTuple, result, mem, config)
    if J === nothing && F !== nothing
        memfNonlinearFunction!(F, x, mee0, meef, tspan, pTuple)
    elseif J !== nothing
        result = ForwardDiff.jacobian!(result, (F,x) -> memfNonlinearFunction!(F,x,mee0,meef,tspan,pTuple), mem, x, config, Val{false}())

        J .= DiffResults.jacobian(result)
        if F !== nothing
            F .= DiffResults.value(result)
        end
    end
    return nothing
end