
using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations
using QLawIndirectOptimization

using Profile
using BenchmarkTools

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
    ephemDays   = 1.0
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
    oe0  = qLawPs.oe0
    mee0 = AstroUtils.convertState(oe0, AstroUtils.Keplerian, AstroUtils.MEE, qLawPs.meePs.mu)

    # QLawIndirectOptimization.qLaw_control(
    #     oe0[1], oe0[2], oe0[3], oe0[4], oe0[5], oe0[6], m0, true, qLawPs 
    # )

    # @code_warntype QLawIndirectOptimization.qLaw_control(
    #     oe0[1], oe0[2], oe0[3], oe0[4], oe0[5], oe0[6], m0, true, qLawPs 
    # )

    @code_warntype QLawIndirectOptimization.qLaw_control_with_coast_checks(mee0, m0, true, qLawPs)

end

main()