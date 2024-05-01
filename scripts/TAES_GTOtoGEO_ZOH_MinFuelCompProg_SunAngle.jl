using AstroEOMs, AstroUtils, SPICE, StaticArrays
using OrdinaryDiffEq
using QLawIndirectOptimization
using Infiltrator

function main()
    furnshDefaults()

    # Compute initial epoch
    initEpoch   = utc2et("2000-03-22T00:00:00")

    # Spacecraft model
    m0   = 1200.0
    P    = 5.0*1000.0   # [W]
    Isp  = 1800.0       # [s]
    g0   = 9.80664      # [m/s^2]
    η    = 0.55 
    tMax = 2*η*P / (g0 * Isp)
    spaceCraft = SimpleSpacecraft(m0, m0, tMax, Isp)

    # Force model parameters
    μs          = 3.986e5
    ephemDays   = 1000.0
    ephemTspan  = (initEpoch - 100.0, initEpoch + ephemDays*86400.0)
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
    mee0        = SA[11359.07, 0.7306, 0.0, 0.2539676, 0.0, 0.0]
    kep0        = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Keplerian, μs)
    kept        = [42165.0, 0.01, 0.01, 0.0, 0.0]

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
        oeW                      = [1.0, 1.0, 1.0, 0.0, 0.0],
        oeTols                   = tolVec,
        ηr_tol                   = 0.1,
        meeParams                = meeParams,
        spaceCraft               = spaceCraft,
        desolver                 = Vern7(),
        maxRevs                  = 500.0,
        integStepOpt             = 10.0,
        integStepGen             = 0.1,
        writeData                = true,
        type                     = :QDSAA,
        eSteps                   = 30,
        eclipsing                = true,
        thrustSunAngleConstraint = true,
        thrustSunAngle           = 50.0*pi/180.0,
        onlyWriteDataAtSteps     = true,
        savedStatesAtSteps       = 10,
    )

    # Define cost
    function cost(state, time, err, retcode)
        J = time - 5*state[7]
        if retcode != :success
            J += 1e8*maximum(err)
        end
        return J
    end

    # Solve
    # cache, meef, kepf, time, retcode = generate_qlaw_transfer(
    #     qLawPs, cost; 
    #     max_time        = 500.0, 
    #     show_trace      = true, 
    #     num_particles   = 100,
    # )
    qLawPs.oeW .= [2.9085806924389734, 2.2829812129445024, 10.0, 0.0, 0.0]
    qLawPs.ηr = 0.18908574785537097
    cache, meef, kepf, time, retcode = generate_qlaw_transfer(qLawPs)
    plot_transfer("test.png", cache, qLawPs)

    @infiltrate
end

main()