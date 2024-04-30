
# Updates the sun location within qLawParams and checks if in eclipse
function update_sun_location!(ps::qLawParams, mee, t)
    ps.eclipsed = false
    if ps.meePs.thirdBodyEphemerides !== nothing && 10 in ps.meePs.thirdBodyEphemerides.targIDs
        # Get position of sun
        r_sun = AstroUtils.getPosition(ps.meePs.thirdBodyEphemerides, 10, ps.meePs.initEpoch + t*ps.meePs.TU)

        # Get position of spacecraft
        cart_us = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Cartesian, ps.μ)
        r_sc = SA[cart_us[1], cart_us[2], cart_us[3]]

        # Get vector from spacecraft to sun
        r_sun_sc = SA[
            r_sun[1] / ps.meePs.LU - r_sc[1],
            r_sun[2] / ps.meePs.LU - r_sc[2],
            r_sun[3] / ps.meePs.LU - r_sc[3],
        ]

        # Compute scale to unit vector and set
        n_r_sun_sc = norm(r_sun_sc)
        ps.toSunVec .= r_sun_sc / n_r_sun_sc

        # Check if in eclipse
        n_r_sc          = norm(r_sc)
        aSR             = asin(ps.RS / n_r_sun_sc)
        aBR             = asin(ps.RB / n_r_sc)
        aD              = acos(-dot(r_sc,r_sun_sc) / (n_r_sc*n_r_sun_sc))
        if ps.eclipsing && aSR + aBR > aD
            ps.eclipsed = true
        end
    end
    return nothing
end

# Computes the error in targeting the desired Keplerian elements
function get_targeting_error(kep, ps)
    aerr = ps.Ws[1]*abs(kep[1] - ps.oet[1]) - ps.oeTols[1]
    eerr = ps.Ws[2]*abs(kep[2] - ps.oet[2]) - ps.oeTols[2]
    ierr = ps.Ws[3]*abs(kep[3] - ps.oet[3]) - ps.oeTols[3]
    Ωerr = ps.Ws[4]*abs(acos(cos(kep[4] - ps.oet[4]))) - ps.oeTols[4]
    ωerr = ps.Ws[5]*abs(acos(cos(kep[5] - ps.oet[5]))) - ps.oeTols[5]
    return (aerr, eerr, ierr, Ωerr, ωerr)
end

# Checks if we've hit any criteria indicating we should stop the sim
function check_stop_criteria(mee, m, error, ps, Lspan)
    if mee[6] >= Lspan[2]
        return :hit_max_L, true
    elseif maximum(error) <= 0.0
        return :success, true
    elseif m < (ps.m0 - ps.mp)
        return :no_fuel, true
    else
        return :continue, false
    end
end

# Performs a single integration step of the qLaw transfer, holding the control constant
# and does not employ a cache (i.e., we don't save any information)
function take_integration_step!(ps, mee, mass, time, L0, Lf, sa_constrained, cache::Nothing)
    # Update sun location
    update_sun_location!(ps, mee, time)
    
    # Compute qLaw control
    α, β, T, e_coast, pq_coast = qLaw_control_with_coast_checks(mee, mass, sa_constrained, ps)

    # Perform numerical integration
    x0 = SA[mee[1], mee[2], mee[3], mee[4], mee[5], time, mass]
    pt = (α, β, T, ps.c, e_coast || pq_coast || ps.eclipsed)
    fn = (u,p,L) -> qLawEOMsSundmanTransformedZOH(u,p,L,ps.meePs)
    sol = solve(
        ODEProblem{false,SciMLBase.FullSpecialize}(fn, x0, (L0,Lf), pt),
        ps.desolver; 
        reltol = ps.reltol,
        abstol = ps.abstol,
        save_everystep  = false,
        save_start      = false,
        initialize_save = false,
        verbose         = false,
    )

    # Form mee state at end of integration
    xf   = sol.u[end]
    meef = SA[xf[1], xf[2], xf[3], xf[4], xf[5], Lf]
    return meef, xf[7], xf[6]
end

# Performs a single integration step of the qLaw transfer, holding the control constant
# and does not employ a cache (i.e., we don't save any information)
function take_integration_step!(ps, mee, mass, time, L0, Lf, sa_constrained, cache::QLawTransferCache)
    # Update sun location
    update_sun_location!(ps, mee, time)
    
    # Compute qLaw control
    α, β, T, eff_coast, pq_coast = qLaw_control_with_coast_checks(mee, mass, sa_constrained, ps)

    # Perform numerical integration
    x0 = SA[mee[1], mee[2], mee[3], mee[4], mee[5], time, mass]
    pt = (α, β, T, ps.c, eff_coast || pq_coast || ps.eclipsed)
    fn = (u,p,L) -> qLawEOMsSundmanTransformedZOH(u,p,L,ps.meePs)
    sol = solve(
        ODEProblem{false,SciMLBase.FullSpecialize}(fn, x0, (L0,Lf), pt),
        ps.desolver; 
        reltol = ps.reltol,
        abstol = ps.abstol,
    )

    # Push update to cache
    push_update!(cache, sol, α, β, T, ps.eclipsed, eff_coast, pq_coast, ps)

    # Form mee state at end of integration
    xf   = sol.u[end]
    meef = SA[xf[1], xf[2], xf[3], xf[4], xf[5], Lf]
    return meef, xf[7], xf[6]
end

function generate_qlaw_transfer(ps::qLawParams, cache::Union{Nothing,QLawTransferCache})
    # Compute initial MEE states
    mee = AstroUtils.convertState(ps.oe0, AstroUtils.Keplerian, AstroUtils.MEE, ps.μ)

    # Initialize time
    time = 0.0

    # Iniitialize mass
    mass = ps.m0

    # Compute independant variable span and info
    Lspan = (mee[6], mee[6] + 2.0*pi*ps.maxRevs)

    # Check if enforcing sa constraint
    sa_con = ps.type != :QDUC

    # Handle step size
    integStep = cache === nothing ? ps.integStepOptimization : ps.integStepGeneration

    # Begin integration loop
    done    = false
    retcode = :continue
    L0      = Lspan[1]
    Lf      = L0 + integStep
    while !done
        mee, mass, time = take_integration_step!(ps, mee, mass, time, L0, Lf, sa_con, cache)

        # Compute kep state at end of integration
        kep = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Keplerian, ps.μ)

        # Compute targeting error
        err = get_targeting_error(kep, ps)

        # Check stopping criteria
        retcode, done = check_stop_criteria(mee, mass, err, ps, Lspan)

        # Update integration span variables
        L0  = Lf 
        Lf  = L0 + integStep
    end

    # Compute final kep state
    kep     = AstroUtils.convertState(mee, AstroUtils.MEE, AstroUtils.Keplerian, ps.μ)
    kepfs   = SA[
        kep[1] * ps.meePs.LU, 
        kep[2], 
        kep[3] * 180.0 / pi,
        kep[4] * 180.0 / pi,  
        kep[5] * 180.0 / pi, 
        kep[6] * 180.0 / pi,
        mass * ps.meePs.MU,
    ]

    # Return final keplerian state, final time, and retcode
    return SA[mee[1],mee[2],mee[3],mee[4],mee[5],Lf,mass], kepfs, time, retcode
end

function generate_qlaw_transfer(ps::qLawParams)
    cache = QLawTransferCache(ps.oeW, ps.ηr)
    meef, kepf, time, retcode = generate_qlaw_transfer(ps, cache)
    return cache, meef, kepf, time, retcode
end

# Function must be of the form:
#   weight_optimization_cost(final_state, final_time, retcode)
#
function generate_qlaw_transfer(
    ps::qLawParams, weight_optimization_cost::F; 
    max_time = 60.0, 
    num_particles = 200, 
    show_trace = false,
) where F <: Function
    # Determine decision variables
    dec_vec_flags = get_decision_variable_flags(ps)
    n = sum(dec_vec_flags)

    # Define optimization problem
    LB = fill(0.0, n)
    UB = fill(10.0, n)
    if dec_vec_flags[end]
        UB[end] = 1.0
    end
    prob = GlobalOptimization.OptimizationProblem(
        x -> qlaw_weight_optimization_cost(x, dec_vec_flags, ps, weight_optimization_cost),
        LB, UB
    )
    pso = GlobalOptimization.PolyesterPSO(
        prob; 
        max_time = max_time, 
        num_particles = num_particles, 
        max_stall_time = Inf, 
        display = show_trace,
    )
    res = GlobalOptimization.optimize!(pso)

    # Update the parameters with the optimal decision variables
    parameter_update_from_dec_vec!(ps, res.xbest, dec_vec_flags)

    # Generate the transfer and return
    return generate_qlaw_transfer(ps)
end

function qlaw_weight_optimization_cost(dec_vec, dec_vec_flags, ps_original, cost::F) where F <: Function
    # Create copy of parameters and update with dec_vec
    ps = deepcopy(ps_original)
    parameter_update_from_dec_vec!(ps, dec_vec, dec_vec_flags)
    
    # Run the simulation
    J = try
        # Run the sim
        mee, kepfs, time, retcode = generate_qlaw_transfer(ps, nothing)

        # Call the user provided cost function
        cost(mee, time, retcode)
    catch e
        if e isa InterruptException
            throw(e)
        else
            1e12
        end
    end
    return J
end

function parameter_update_from_dec_vec!(ps, dec_vec, dec_vec_flags)
    idx = 1
    for i in eachindex(ps.oeW)
        if dec_vec_flags[i]
            ps.oeW[i] = dec_vec[idx]
            idx += 1
        end
    end
    if dec_vec_flags[6]
        ps.ηr = dec_vec[end]
    end
end

function get_decision_variable_flags(ps::qLawParams)
    # Define decision variable flags
    flag_a   = ps.oeW[1] > 0.0
    flag_e   = ps.oeW[2] > 0.0
    flag_i   = ps.oeW[3] > 0.0
    flag_ω   = ps.oeW[4] > 0.0
    flag_Ω   = ps.oeW[5] > 0.0
    flag_eff = ps.ηr > 0.0
    return (flag_a, flag_e, flag_i, flag_ω, flag_Ω, flag_eff)
end

