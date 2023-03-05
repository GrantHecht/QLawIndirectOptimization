# meef_interp is a container of interpolants for geting final state constraint and time
function minFuelMayerPseudoarclengthNonlinearFunction!(F, x, xi, dc, Δs, mee0, meef_interp, pTuple)
    # Compute optimal control nonlinear function
    @views minFuelMayerNonlinearFunctionWithInterp!(F[1:6], x, mee0, meef_interp, pTuple)

    # Compute step constraint
    F[7] = -Δs
    @inbounds for i in 1:7
        F[7] += (x[i] - xi[i])*dc[i]
    end
end

function minFuelMayerPseudoarclengthNonlinearFunctionWithJacobian!(F, J, x, xi, dc, Δs, mee0, meef_interp, pTuple, result, mem, config)
    if J === nothing && F !== nothing
        minFuelMayerPseudoarclengthNonlinearFunction!(F, x, xi, dc, Δs, mee0, meef_interp, pTuple)
    elseif J !== nothing
        result = ForwardDiff.jacobian!(result, 
                    (F,x) -> minFuelMayerPseudoarclengthNonlinearFunction!(F, x, xi, dc, Δs, mee0, meef_interp, pTuple), 
                    mem, x, config, Val{false}())

        J .= DiffResults.jacobian(result)
        if F !== nothing
            F .= DiffResults.value(result)
        end
    end
    return nothing
end

# meef_interp is a container of interpolants for geting final state constraint and time
function minFuelMayerNonlinearFunctionWithInterp!(F, x, mee0, meef_interp, pTuple)
    # Grab the interpolated state for the current value of κ 
    val = getFinalStateAndTime(meef_interp, x[7])

    # Construct full initial state
    y0  = SVector(mee0[1], mee0[2], mee0[3], mee0[4], mee0[5], mee0[6], pTuple[1].initMass / pTuple[2].MU,
                    x[1], x[2], x[3], x[4], x[5], x[6], -1.0)

    # Construct problem and solve
    prob = ODEProblem(ODEFunction{false}(AstroEOMs.meeOptimalControlEOMs), y0, (0.0, val[6]), pTuple)
    sol  = solve(prob, Vern9(), reltol = 1e-14, abstol = 1e-14, save_everystep = false, save_start = false, initialize_save = false, maxiters = 1e7,
            sensealg = ForwardDiffSensitivity(convert_tspan=true))
    xf   = sol[end]

    F[1] = xf[1] - val[1]
    F[2] = xf[2] - val[2] 
    F[3] = xf[3] - val[3]
    F[4] = xf[4] - val[4]
    F[5] = xf[5] - val[5]
    F[6] = xf[13]
end
