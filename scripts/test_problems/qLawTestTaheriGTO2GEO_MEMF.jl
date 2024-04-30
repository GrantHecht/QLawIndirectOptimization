using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
furnshDefaults()

function main()
    # Compute initial epoch
    initEpoch   = utc2et("2012-12-23T12:00:00")

    # Define parameters for EOMs
    μs          = 3.986e5
    tMax        = 0.5
    Isp         = 3100.0
    m0          = 100.0
    mp          = 90.0
    LU          = 6378.0
    TU          = sqrt(LU^3 / μs)
    meeParams   = MEEParams(initEpoch; LU = LU, MU = 1.0, TU = TU, μ = μs, 
                    costFunction = AstroEOMs.MinimumEnergyToFuel)
    spaceCraft  = SimpleSpacecraft(m0, mp, tMax, Isp)
    sw          = SwitchingStruct(1.0, 1)

    # Define initial and target orbital elements
    cart0       = SVector(6378.9, 0.0, 0.0,
                    0.0, 10.0258, 1.231)
    kep0,f      = AstroUtils.convertState(cart0, AstroUtils.Cartesian, AstroUtils.Keplerian, μs)
    mee0        = Vector(AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μs))
    kep0d       = Vector(kep0)
    for i in 3:6
        kep0d[i] = mod(kep0d[i], 2.0*pi) * 180.0 / pi
    end
    kept        = [42165.0, 0.01, 0.01, 0.0, 0.0]

    # Define qLaw parameters
    oeW          = [1.0, 1.0, 5.0, 0.0, 0.0] 

    # Define tolerance on targeted elements
    atol        = 20.0
    etol        = 0.001
    itol        = 0.01
    Ωtol        = 0.01
    ωtol        = 0.01
    tolVec      = [atol,etol,itol,Ωtol,ωtol]

    # Construct qLaw parameters
    qLawPs       = qLawParams(copy(kep0d), copy(kept);
                    oeW         = oeW,
                    oeTols      = tolVec,
                    ηr_tol      = 0.15,
                    ηa_tol      = 0.0,
                    meeParams   = meeParams,
                    spaceCraft  = spaceCraft,
                    desolver    = Vern7(),
                    maxRevs     = 20.0,
                    integStep   = 2.0*pi / 180.0,
                    writeData   = true,
                    returnData  = true,
                    type        = :QDUC)

    # Run QLaw sim
    meeAtSteps, kepf, retcode = qLawOriginal(qLawPs)

    # Scale initial mee state
    mee0[1] /= meeParams.LU

    # Grab state at ΔL = 2*pi
    L0      = meeAtSteps[1,6]
    istart  = 0
    for i in axes(meeAtSteps,1)
        if abs(meeAtSteps[i,6] - L0 - 2*pi) < 1e-3
            istart = i
            break
        end
    end
    meef     = meeAtSteps[istart,1:7]
    tf       = meeAtSteps[istart,8]

    # Solve short problem
    λ0, retcode = memfSolve(mee0, meef, (0.0, tf), (spaceCraft,meeParams,sw))

    # Compute percent of trajectory which we've solved for
    percs  = floor(Int,100.0 * (istart / size(meeAtSteps, 1)))
    println(string(percs) * "% complete.")

    # If successful, continue continuation along qLaw trajectory
    ftol_steps = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6]
    λ0g     = zeros(7)
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
                λ0, retcode = memfSolve(mee0, meef, (0.0, tf), (spaceCraft,meeParams,sw); 
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
                idx += 10
            else
                # Solve again and print trace
                λ0, retcode = memfSolve(mee0, meef, (0.0, tf), (spaceCraft,meeParams,sw); 
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
        open(datadir("TaheriGTO2GEO_MEMF_costates_new.txt"), "w") do io; writedlm(io, λ0sm); end
        open(datadir("TaheriGTO2GEO_MEMF_finalStates_new.txt"), "w") do io; writedlm(io, meefsm); end
        open(datadir("TaheriGTO2GEO_MEMF_tofs_new.txt"), "w") do io; writedlm(io, tfs); end
        open(datadir("TaheriGTO2GEO_MEMF_ftols_new.txt"), "w") do io; writedlm(io, ftols); end
        open(datadir("TaheriGTO2GEO_MEMF_percs_new.txt"), "w") do io; writedlm(io, percents); end
    catch 
        println("Failed to write files.")
    end

    return (λ0s, meefs, tfs, ftols, percents)
end

λ0s, meefs, tfs, ftols, percents = main()

