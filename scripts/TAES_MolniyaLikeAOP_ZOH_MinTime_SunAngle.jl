using AstroEOMs, AstroUtils, SPICE, StaticArrays
using OrdinaryDiffEq
using QLawIndirectOptimization
using JLD2, Infiltrator

function main()
    # Furnsh default spice kernels
    furnshDefaults()

    # Compute initial epoch
    initEpoch   = utc2et("2023-01-01T13:30:00")

    # Spacecraft model
    m0   = 1200.0
    Isp  = 3000.0    
    tMax = 0.5
    spaceCraft = SimpleSpacecraft(m0, m0, tMax, Isp)

    # Force model parameters
    μs          = 3.986e5
    ephemDays   = 5000.0
    nPoints     = ceil(Int64, 2*ephemDays)
    tbEphems    = Ephemerides(
        (initEpoch - 100.0, initEpoch + ephemDays*86400.0), 
        nPoints, 
        [10], 
        399, 
        "J2000",
    )
    meeParams   = MEEParams(
        initEpoch; 
        LU = 1.0, MU = 1.0, TU = 24.0*3600.0, 
        μ = μs,
        thirdBodyPerterbations = true,
        thirdBodyEphemerides = tbEphems, 
        nonsphericalGravity = true,
    )

    # Define initial and target orbital elements
    mee0        = SA[11359.07, 0.7306, 0.0, 0.2539676, 0.0, 0.0]
    kep0        = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Keplerian, μs)
    kept        = [26500.0, 0.700, 116.0, 180.0, 270.0]

    # Convert angles in initial kep state to deg
    kep0d       = [kep0[i] for i in eachindex(kep0)]
    kep0d[3:6] .*= 180.0 / pi

    # Define tolerance on targeted elements
    atol        = 20.0
    etol        = 0.001
    itol        = 0.01
    Ωtol        = 0.01
    ωtol        = 0.01
    tolVec      = [atol,etol,itol,Ωtol,ωtol]

    # Construct qLaw parameters
    qLawPs       = qLawParams(
        kep0d, kept;
        oeW                      = [1.0, 1.0, 1.0, 0.0, 1.0],
        Wp                       = 0.1,
        rpmin                    = 6578.0,
        oeTols                   = tolVec,
        ηr_tol                   = 0.0,
        meeParams                = meeParams,
        spaceCraft               = spaceCraft,
        desolver                 = Vern7(),
        reltol                   = 1e-8,
        abstol                   = 1e-8,
        maxRevs                  = 600.0,
        integStepOpt             = 1.0,
        integStepGen             = 1.0,
        writeData                = true,
        type                     = :QDSA,
        eSteps                   = 10,
        eclipsing                = true,
        thrustSunAngleConstraint = true,
        thrustSunAngle           = 50.0*pi/180.0,
        onlyWriteDataAtSteps     = true,
        savedStatesAtSteps       = 10,
    )

    # Define cost
    function cost(state, time, err, retcode)
        J = time #- 5*state[7]
        if retcode != :success
            J = 1e8*maximum(err)
        end
        return J
    end

    # Solve
    cache, meef, kepf, time, retcode = generate_qlaw_transfer(
        qLawPs, cost, QLawIndirectOptimization.ThreadedPSO; 
        max_time        = 12*3600.0, 
        show_trace      = true, 
        num_particles   = 200,
    )

    # From 12 hour optimization
    qLawPs.oeW .= [10.0,1.2089577606081345,0.873052700670224,0.0,0.01633465971623084] 
    cache, meef, kepf, time, retcode = generate_qlaw_transfer(qLawPs)

    # qLawPs.oeW .= [9.546535123066946, 4.9366795064551185, 3.360605913416598, 1.202994187072187, 0.0]
    # cache, meef, kepf, time, retcode = generate_qlaw_transfer(qLawPs)

    # Save solution information
    jldsave(
        joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLikeAOP_MinTime.jld2");
        cache = cache, params = qLawPs,
    )
    dump_to_mat(cache, joinpath(@__DIR__, "..", "data", "TAES", "mat", "MolniyaLikeAOP_MinTime.mat"))

    # Generate figures
    plot_transfer(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "MolniyaLikeAOP_MinTime_xy.png"), cache, qLawPs; 
        axes = [1,2], linewidth=0.2,
    )
    plot_transfer(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "MolniyaLikeAOP_MinTime_xz.png"), cache, qLawPs; 
        axes = [1,3], linewidth=0.2,
    )
    plot_transfer(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "MolniyaLikeAOP_MinTime_yz.png"), cache, qLawPs; 
        axes = [2,3], linewidth=0.2,
    )

    @infiltrate
    return (cache, qLawPs)
end

cache, ps = main()