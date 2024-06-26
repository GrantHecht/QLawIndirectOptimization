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
    tMax = 2.0
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
        thirdBodyPerterbations = false,
        thirdBodyEphemerides = tbEphems, 
        nonsphericalGravity = false,
    )

    # Define initial and target orbital elements
    kep0        = [24505.9, 0.725,   0.6,   0.0,   0.0, 0.0]
    kept        = [26500.0, 0.700, 116.0, 180.0, 270.0]

    # Define tolerance on targeted elements
    atol        = 20.0
    etol        = 0.001
    itol        = 0.01
    Ωtol        = 0.01
    ωtol        = 0.01
    tolVec      = [atol,etol,itol,Ωtol,ωtol]

    # Construct qLaw parameters
    qLawPs       = qLawParams(
        kep0, kept;
        oeW                      = [1.0, 1.0, 1.0, 1.0, 1.0],
        Wp                       = 1.0,
        rpmin                    = 6578.0,
        oeTols                   = tolVec,
        ηr_tol                   = 0.0,
        meeParams                = meeParams,
        spaceCraft               = spaceCraft,
        desolver                 = Tsit5(),
        reltol                   = 1e-8,
        abstol                   = 1e-8,
        maxRevs                  = 200.0,
        integStepOpt             = 10.0,
        integStepGen             = 0.1,
        writeData                = true,
        type                     = :QDUC,
        eSteps                   = 10,
        eclipsing                = true,
        thrustSunAngleConstraint = true,
        thrustSunAngle           = 50.0*pi/180.0,
        onlyWriteDataAtSteps     = true,
        savedStatesAtSteps       = 10,
    )

    # Define cost
    function cost(state, time, retcode)
        J = time #- 5*state[7]
        if retcode != :success
            J += 1e12
        end
        return J
    end

    # Solve
    cache, meef, kepf, time, retcode = generate_qlaw_transfer(
        qLawPs, cost, QLawIndirectOptimization.PolyesterPSO; 
        max_time        = 3600.0, 
        show_trace      = true, 
        num_particles   = 200,
    )
    # cache, meef, kepf, time, retcode = generate_qlaw_transfer(qLawPs)
    plot_transfer("test.png", cache, qLawPs; axes = SA[1,2])

    # Save
    jldsave(
        joinpath(@__DIR__, "..", "data", "TAES", "PetroMolniya.jld2"), 
        cache = cache, qLawPs = qLawPs,
    )

    @infiltrate
    return (cache, qLawPs)
end

cache, ps = main()