using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
using Heuristics
furnshDefaults()

function main()
    # Compute initial epoch
    initEpoch   = utc2et("2022-10-07T12:00:00")

    # Compute thrust used by Shannon et al.
    P    = 5.0*1000.0   # [W]
    Isp  = 1800.0       # [s]
    g0   = 9.80664      # [m/s^2]
    η    = 0.55 
    tMax = 2*η*P / (g0 * Isp)

    # Define parameters for EOMs
    μs          = 3.986e5
    m0          = 1200.0
    LU          = 6378.0
    TU          = sqrt(LU^3 / μs)
    meeParams   = MEEParams(initEpoch; LU = LU, MU = 1.0, TU = TU, μ = μs,
                    costFunction = AstroEOMs.MinimumFuelMayer)
    spaceCraft  = SimpleSpacecraft(m0, m0, tMax, Isp)
    sw          = SwitchingStruct(1.0, 1)

    # Define initial and target orbital elements
    mee0        = [11359.07, 0.7306, 0.0, 0.2539676, 0.0, 0.0]
    kep0, f     = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Keplerian, μs)
    kept        = [42165.0, 0.01, 0.01, 0.0, 0.0]

    # Convert angles in initial kep state to deg
    kep0d        = Vector(kep0)
    kep0d[3:6] .*= 180.0 / pi

    # Define error weights
    oeW         = [1.193, 2.402, 8.999, 0.0, 0.0] 

    # Define tolerance on targeted elements
    atol        = 20.0
    etol        = 0.001
    itol        = 0.01
    Ωtol        = 0.01
    ωtol        = 0.01
    tolVec      = [atol,etol,itol,Ωtol,ωtol]

    # Construct qLaw parameters
    qLawPs       = qLawParams(kep0d, kept;
                    oeW         = oeW,
                    oeTols      = tolVec,
                    ηr_tol      = 0.1,
                    meeParams   = meeParams,
                    spaceCraft  = spaceCraft,
                    desolver    = Vern7(),
                    maxRevs     = 1000.0,
                    integStep   = 5.0,
                    writeData   = true,
                    returnData  = true,
                    type        = :QDUC)

    # Run QLaw sim
    meeAtSteps, kepf, retcode = qLawOriginal(qLawPs, :minfuel)

    # Scale initial mee state
    mee0[1] /= meeParams.LU

    # Grab state at ΔL = 2*pi
    L0      = meeAtSteps[1,6]
    istart  = 0
    for i in axes(meeAtSteps,1)
        if abs(meeAtSteps[i,6] - L0 - 6.0*pi) < 1e-3
            istart = i
            break
        end
    end
    meef     = meeAtSteps[istart,1:7]
    tf       = meeAtSteps[istart,8]

    # Solve short problem
    λ0, retcode = minFuelMayerSolve(mee0, meef, (0.0, tf), (spaceCraft,meeParams,sw); 
                    show_trace = true, ftol = 1e-8)

    # Compute percent of trajectory which we've solved for
    percs  = floor(Int,100.0 * (istart / size(meeAtSteps, 1)))
    println(string(percs) * "% complete.")

    # If successful, continue continuation along qLaw trajectory
    ftol_steps = [5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6]
    λ0g     = zeros(6)
    λ0s     = [λ0]
    meefs   = [meef]
    tfs     = [tf]
    lastFtol= 0.0
    ftols   = [1e-8]
    percents= [100.0 * (istart / size(meeAtSteps, 1))]
    if retcode == :success
        idx = istart + 1
        while idx <= size(meeAtSteps, 1)
            # Print out percent complete if we've improved by a while pecent
            perc  = 100.0 * (idx / size(meeAtSteps, 1))
            push!(percents, perc)
            if perc >= percs + 1.0
                percs = floor(Int, perc)
                println(string(percs) * "% complete.")
            end

            # Grab new final state and time
            meef .= meeAtSteps[idx, 1:7]
            tf = meeAtSteps[idx, 8]

            # Solve using previous solution as guess with decreasing ftol
            success = true
            λ0g    .= λ0s[end]
            run(`clear`)
            println(string(percs) * "% complete.")
            for ftol in ftol_steps
                print("Solving with ftol = " * string(ftol) * "...")
                λ0, retcode = minFuelMayerSolve(mee0, meef, (0.0, tf), (spaceCraft,meeParams,sw); 
                                λ0g = λ0g, ftol = ftol, iterations = 200)
                if retcode == :success
                    print(" Success!\n")
                    λ0g .= λ0
                    lastFtol = ftol
                else
                    print(" Failed.\n")
                    if ftol == ftol_steps[1]
                        success = false
                    end
                    break
                end
            end

            # Update idx if successful
            if success
                push!(λ0s, copy(λ0g))
                push!(meefs, copy(meef))
                push!(tfs, tf)
                push!(ftols, lastFtol)
                idx += 1
            else
                # Solve again and print trace
                λ0, retcode = minFuelMayerSolve(mee0, meef, (0.0, tf), (spaceCraft,meeParams,sw); 
                                λ0g = λ0s[end], show_trace = true)
                break
            end
        end
    end

    # Write return values to file
    try
        λ0sm    = zeros(length(λ0s), length(λ0s[1]))
        for i in eachindex(λ0s)
            λ0sm[i,:] .= λ0s[i]
        end
        meefsm  = zeros(length(meefs), length(meefs[1]))
        for i in eachindex(meefs)
            meefsm[i,:] .= meefs[i]
        end
        open(datadir("Shannon_costates_new.txt"), "w") do io; writedlm(io, λ0sm); end
        open(datadir("Shannon_finalStates_new.txt"), "w") do io; writedlm(io, meefsm); end
        open(datadir("Shannon_tofs_new.txt"), "w") do io; writedlm(io, tfs); end
        open(datadir("Shannon_ftols_new.txt"), "w") do io; writedlm(io, ftols); end
        open(datadir("Shannon_percs_new.txt"), "w") do io; writedlm(io, percents); end
    catch 
        println("Failed to write files.")
    end

    return (λ0s, meefs, tfs, ftols, percents)
end

λ0s, meefs, tfs, ftols, percents = main()
