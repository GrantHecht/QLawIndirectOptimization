using AstroEOMs, AstroUtils, SPICE, StaticArrays
using OrdinaryDiffEq
using QLawIndirectOptimization
using Infiltrator

function main()
    # Furnsh default spice kernels
    furnshDefaults()

    # Compute initial epoch
    initEpoch   = utc2et("2023-01-01T13:30:00")

    # Spacecraft model
    m0   = 1200.0
    Isp  = 1800.0    
    tMax = 0.3
    spaceCraft = SimpleSpacecraft(m0, m0, tMax, Isp)

    # Force model parameters
    μs          = 3.986e5
    ephemDays   = 5000.0
    nPoints     = ceil(Int64, 2*ephemDays)
    tbEphems    = Ephemerides(
        (initEpoch - 100.0, initEpoch + ephemDays*86400.0), 
        nPoints, 
        [10,301], 
        399, 
        "J2000",
    )
    meeParams   = MEEParams(
        initEpoch; 
        LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, 
        μ = μs,
        thirdBodyEphemerides = tbEphems, 
        nonsphericalGravity = true,
    )

    # Define initial and target orbital elements
    kep0        = [26500.0, 0.7, 0.01,  0.0,  0.0, 0.0]
    kept        = [26700.0, 0.7, 20.0, 60.0,  0.0]

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
        oeW                      = [1.0, 1.0, 1.0, 1.0, 0.0],
        oeTols                   = tolVec,
        ηr_tol                   = 0.0,
        meeParams                = meeParams,
        spaceCraft               = spaceCraft,
        desolver                 = Vern7(),
        maxRevs                  = 500.0,
        integStepOpt             = 10.0,
        integStepGen             = 1.0,
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
    function cost(state, time, err, retcode)
        J = time #- 5*state[7]
        if retcode != :success
            J += 1e8*maximum(err)
        end
        return J
    end

    # Solve
    # cache, meef, kepf, time, retcode = generate_qlaw_transfer(
    #     qLawPs, cost, QLawIndirectOptimization.ThreadedPSO; 
    #     max_time        = 3600.0, 
    #     show_trace      = true, 
    #     num_particles   = 100,
    # )
    qLawPs.oeW .= [4.252647196097554,7.828168230326815,5.644962073810314,0.20125849969261253,0.0]
    cache, meef, kepf, time, retcode = generate_qlaw_transfer(qLawPs)
    plot_transfer("short_xy.png", cache, qLawPs; axes = SA[1,2], linewidth = 0.25)
    plot_transfer("short_xz.png", cache, qLawPs; axes = SA[1,3], linewidth = 0.25)
    plot_transfer("short_yz.png", cache, qLawPs; axes = SA[2,3], linewidth = 0.25)

    @infiltrate
    return (cache, qLawPs)
end

cache, ps = main()