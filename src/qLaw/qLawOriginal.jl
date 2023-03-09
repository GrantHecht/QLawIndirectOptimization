function qLawOriginalCost(x, psOG, decVecFlags, cost)
    # Create copy to remain thread safe
    ps = deepcopy(psOG)

    # Set decision variables in parameters
    idx = 1
    for i in eachindex(ps.oeW)
        if decVecFlags[i] == true
            ps.oeW[i] = x[idx]
            idx += 1
        end
    end
    if decVecFlags[7] == true
        ps.ηr = x[end]
    end

    # Run simulation
    J = try
        meefs, kepfs, retcode = qLawOriginal(ps)

        # Compute cost
        J = 0.0
        if cost == :minfuel
            J -= meefs[7]
        end
        if retcode != :success
            J += 5000.0
        end
        J
    catch 
        10000.0
    end
    return J
end

function qLawOriginal(ps::qLawParams, cost; numParticles = 50, useParallel = true)
    # Get values of returnTrajAtSteps and writeDataToFile before changing
    returnDataOG    = ps.returnTrajAtSteps
    writeDataOG     = ps.writeDataToFile
    ps.returnTrajAtSteps = false
    ps.writeDataToFile   = false

    # Get the number of parameters we're optimizing
    decVecFlags     = [false, false, false, false, false, false, false, false]
    for i in eachindex(ps.oeW)
        if ps.oeW != 0.0
            decVecFlags[i] = true
        end
    end
    if ps.ηr != 0.0
        decVecFlags[7] = true
    end

    # Set decision vector bounds (for full dec vec; will be reduced if not targeting
    # specific elements or using effectivity)
    decVecUB        = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.8]
    decVecLB        = zeros(7)

    # Get decision vector bounds for current problem
    N               = sum(decVecFlags)
    decVecLBC       = zeros(N)
    decVecUBC       = zeros(N)
    idx             = 1
    for i in eachindex(decVecFlags)
        if decVecFlags[i] == true
            decVecUBC[idx] = decVecUB[i]
            idx += 1
        end
    end

    # Optimize
    prob    = Heuristics.Problem(x -> qLawOriginalCost(x, ps, decVecFlags, cost), decVecLBC, decVecUBC)
    opts    = Heuristics.Options(;display = true, useParallel = useParallel)
    pso     = Heuristics.PSO(prob; numParticles = numParticles)
    Heuristics.optimize!(pso, opts)

    # Reset parameters, simulate, and return
    idx = 1
    for i in eachindex(ps.oeW)
        if decVecFlags[i] == true
            ps.oeW[i] = pso.results.xbest[idx]
            idx += 1
        end
    end
    if decVecFlags[7] == true
        ps.ηr = pso.results.xbest[end]
    end
    ps.returnTrajAtSteps = returnDataOG
    ps.writeDataToFile   = writeDataOG

    return qLawOriginal(ps)
end

function qLawOriginal(ps::qLawParams)
    # Compute initial MEE states
    mee0 = AstroUtils.convertState(ps.oe0, AstroUtils.Keplerian, 
                AstroUtils.MEE, ps.μ)

    # Construct initial sundman transformed state
    x0  = SVector(mee0[1], mee0[2], mee0[3], mee0[4], mee0[5], 0.0, ps.m0)

    # Compute independant variable span and info
    Lspan = (mee0[6], mee0[6] + 2.0*pi*ps.maxRevs)
    
    # Allocate storage if saving data
    if ps.writeDataToFile && ps.writeDataOnlyAtSteps == false
        n           = ceil(Int, 360.0*ps.maxRevs)
        Ls          = range(Lspan[1], Lspan[2]; length = n)
        ts          = fill(NaN, n)
        cart_th     = fill(NaN, n, 7)
        mee_th      = fill(NaN, n, 7)
        kep_th      = fill(NaN, n, 7)
        coast_th    = fill(true, n)
        angles_th   = fill(NaN, n, 2)
        thrust_th   = fill(NaN, n)
    elseif ps.writeDataToFile
        ts              = Vector{Float64}(undef, 0)
        mee_th          = Vector{Vector{Float64}}(undef, 0)
        cart_th         = Vector{Vector{Float64}}(undef, 0)
        kep_th          = Vector{Vector{Float64}}(undef, 0)
        coast_th        = Vector{Bool}(undef, 0)
        angles_th       = Vector{Vector{Float64}}(undef, 0)
        thrust_th       = Vector{Float64}(undef, 0)
    end
    if ps.returnTrajAtSteps
        ns                  = ceil(Int, (2.0*pi / ps.integStep)*ps.maxRevs)
        mee_th_steps        = fill(NaN, ns, 8)
        mee_th_steps[1,1]   = mee0[1]
        mee_th_steps[1,2]   = mee0[2]
        mee_th_steps[1,3]   = mee0[3]
        mee_th_steps[1,4]   = mee0[4]
        mee_th_steps[1,5]   = mee0[5]
        mee_th_steps[1,6]   = mee0[6]
        mee_th_steps[1,7]   = ps.m0
        mee_th_steps[1,8]   = 0.0
    end

    # Begin integration loop
    done    = false
    idx     = 1
    idxs    = 1
    L0      = Lspan[1]
    Lf      = L0 + ps.integStep
    retcode = :none
    while !done
        # Construct MEE state
        mee     = SVector(x0[1], x0[2], x0[3], x0[4], x0[5], L0)

        # Compute keplerian state
        kep,f   = AstroUtils.convertState(mee, AstroUtils.MEE,
                    AstroUtils.Keplerian, ps.μ)
        if f != 0
            retcode = :keplerian_singularity
            break
        end

        # Compute qLaw control
        α,β,T,coast = qLawThrust_Keplerian(mee, x0[7], ps; method = ps.type)
        ps.α  = α
        ps.β  = β
        ps.T  = T
        ps.coasting = coast

        # Perform numerical integration
        prob    = ODEProblem(qLawEOMsSundmanTransformedZOH, x0, (L0,Lf), ps)
        sol     = solve(prob, ps.desolver, reltol = ps.reltol, abstol = ps.abstol)

        # Save info if desired
        if ps.writeDataToFile && ps.writeDataOnlyAtSteps == false
            while idx <= n && Ls[idx] <= sol.t[end] && Ls[idx] >= sol.t[1] 
                mees_us     = sol(Ls[idx])
                mee_us      = SVector(mees_us[1], mees_us[2], mees_us[3],
                                mees_us[4], mees_us[5], Ls[idx], mees_us[7])

                ts[idx]             = mees_us[6]
                mee_th[idx,1]       = ps.meePs.LU * mee_us[1]
                mee_th[idx,2:6]    .= mee_us[2:6] 
                mee_th[idx,7]       = ps.meePs.MU * mee_us[7]

                cart_us             = AstroUtils.convertState(mee_us, AstroUtils.MEE, 
                                        AstroUtils.Cartesian, ps.μ)
                cart_th[idx,1:3]   .= ps.meePs.LU*view(cart_us,1:3)
                cart_th[idx,4:6]   .= ps.meePs.LU*view(cart_us,4:6)/ps.meePs.TU
                cart_th[idx,7]      = ps.meePs.MU * mee_us[7]

                kep_us,fff          = AstroUtils.convertState(cart_us, AstroUtils.Cartesian, 
                                        AstroUtils.Keplerian, ps.μ)
                kep_th[idx,1]       = ps.meePs.LU*kep_us[1]
                kep_th[idx,2:6]    .= view(kep_us, 2:6)
                kep_th[idx,7]       = ps.meePs.MU * mee_us[7]

                # Deal with coasting
                coast_th[idx]      = ps.coasting

                # Deal with angles
                angles_th[idx,1]   = ps.α
                angles_th[idx,2]   = ps.β
                thrust_th[idx]     = ps.coasting ? 0.0 : ps.T

                # Increment index
                idx += 1
            end
        elseif ps.writeDataToFile
            mees_us     = sol[1]
            mee_us      = SVector(mees_us[1], mees_us[2], mees_us[3],
                            mees_us[4], mees_us[5], Lf, mees_us[7])

            push!(ts, mees_us[6])
            push!(mee_th, zeros(7))
            mee_th[end][1]      = ps.meePs.LU * mee_us[1]
            mee_th[end][2:6]   .= mee_us[2:6] 
            mee_th[end][7]      = ps.meePs.MU * mee_us[7]

            push!(cart_th, zeros(7))
            cart_us             = AstroUtils.convertState(mee_us, AstroUtils.MEE, 
                                    AstroUtils.Cartesian, ps.μ)
            cart_th[end][1:3]  .= ps.meePs.LU*view(cart_us,1:3)
            cart_th[end][4:6]  .= ps.meePs.LU*view(cart_us,4:6)/ps.meePs.TU
            cart_th[end][7]     = ps.meePs.MU * mee_us[7]

            push!(kep_th, zeros(7))
            kep_us,fff          = AstroUtils.convertState(cart_us, AstroUtils.Cartesian, 
                                    AstroUtils.Keplerian, ps.μ)
            kep_th[end][1]      = ps.meePs.LU*kep_us[1]
            kep_th[end][2:6]   .= view(kep_us, 2:6)
            kep_th[end][7]      = ps.meePs.MU * mee_us[7]

            # Deal with coasting
            push!(coast_th, ps.coasting)

            # Deal with angles
            push!(angles_th, [ps.α,ps.β])
            push!(thrust_th, ps.coasting ? 0.0 : ps.T)
        end

        # Save data at step if desired
        if ps.returnTrajAtSteps
            mee_th_steps[idxs,1] = sol[end][1]
            mee_th_steps[idxs,2] = sol[end][2]
            mee_th_steps[idxs,3] = sol[end][3]
            mee_th_steps[idxs,4] = sol[end][4]
            mee_th_steps[idxs,5] = sol[end][5]
            mee_th_steps[idxs,6] = Lf
            mee_th_steps[idxs,7] = sol[end][7]
            mee_th_steps[idxs,8] = sol[end][6]
            idxs += 1
        end
        
        # Compute targeting error
        aerr                = ps.Ws[1]*abs(kep[1] - ps.oet[1]) - ps.oeTols[1]
        eerr                = ps.Ws[2]*abs(kep[2] - ps.oet[2]) - ps.oeTols[2]
        ierr                = ps.Ws[3]*abs(kep[3] - ps.oet[3]) - ps.oeTols[3]
        Ωerr                = ps.Ws[4]*abs(acos(cos(kep[4] - ps.oet[4]))) - ps.oeTols[4]
        ωerr                = ps.Ws[5]*abs(acos(cos(kep[5] - ps.oet[5]))) - ps.oeTols[5]
        targError           = (aerr, eerr, ierr, Ωerr, ωerr)

        # Check stopping criteria
        if Lf >= Lspan[2]
            retcode = :hit_max_L
            done = true
        elseif maximum(targError) <= 0.0
            retcode = :success
            done = true
        elseif sol[end][7] < (ps.m0 - ps.mp)
            retcode = :no_fuel
            done = true
        end

        # Update loop variables
        x0  = SVector(sol[end][1], sol[end][2], sol[end][3], sol[end][4],
                sol[end][5], sol[end][6], sol[end][7])
        L0  = Lf 
        Lf  = L0 + ps.integStep
    end

    # Write to file if desired
    if ps.writeDataToFile
        # Construct target orbit states with SI units
        kepts           = [ps.meePs.LU*ps.oet[1], ps.oet[2], ps.oet[3], 
                            ps.oet[4], ps.oet[5]]

        # Construct constants output
        consts          = [ps.μ * ps.meePs.LU^3 / ps.meePs.TU^2]

        # Write data to files
        open(datadir("kep.txt"),   "w") do io; writedlm(io,   kep_th); end
        open(datadir("mee.txt"),   "w") do io; writedlm(io,   mee_th); end
        open(datadir("cart.txt"),  "w") do io; writedlm(io,  cart_th); end
        open(datadir("coast.txt"), "w") do io; writedlm(io, Int.(coast_th)); end
        open(datadir("angles.txt"),"w") do io; writedlm(io, angles_th); end
        open(datadir("thrust.txt"),"w") do io; writedlm(io, thrust_th); end
        open(datadir("time.txt"),  "w") do io; writedlm(io, ts); end
        open(datadir("kept.txt"),  "w") do io; writedlm(io, kepts); end
        open(datadir("consts.txt"),"w") do io; writedlm(io, consts); end
    end

    # Compute final kep state
    mee     = SVector(x0[1], x0[2], x0[3], x0[4], x0[5], Lf)
    kepf,f  = AstroUtils.convertState(mee, AstroUtils.MEE,
                AstroUtils.Keplerian, ps.μ)
    kepfs   = SVector(kepf[1] * ps.meePs.LU, kepf[2], kepf[3] * 180.0 / pi,
                kepf[4] * 180.0 / pi, kepf[5] * 180.0 / pi, kepf[6] * 180.0 / pi,
                x0[7] * ps.meePs.MU)

    # If returning state at steps, cleanup array before returning
    if ps.returnTrajAtSteps
        stopidx = 0
        for i in axes(mee_th_steps,1)
            if isnan(mee_th_steps[i,1])
                stopidx = i 
                break
            end
        end
        mee_th_steps = mee_th_steps[1:stopidx - 1, :]
    end

    # Return final keplerian state, final time, and retcode
    if ps.returnTrajAtSteps
        return (mee_th_steps, kepfs, retcode)
    else
        return ([x0[1] x0[2] x0[3] x0[4] x0[5] Lf x0[7]], kepfs, retcode)
    end
end