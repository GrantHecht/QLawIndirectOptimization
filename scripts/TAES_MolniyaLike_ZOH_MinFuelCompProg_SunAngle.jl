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
        [10,301], 
        399, 
        "J2000",
    )
    meeParams   = MEEParams(
        initEpoch; 
        LU = 1.0, MU = 1.0, TU = 24.0*3600.0, 
        μ = μs,
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
        oeW                      = [1.0, 1.0, 1.0, 1.0, 0.0],
        Wp                       = 0.1,
        rpmin                    = 6578.0,
        oeTols                   = tolVec,
        ηr_tol                   = 0.1,
        meeParams                = meeParams,
        spaceCraft               = spaceCraft,
        desolver                 = Vern7(),
        reltol                   = 1e-10,
        abstol                   = 1e-10,
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
        J = time - 5*state[7]
        if retcode != :success
            J = 1e8*maximum(err)
        end
        return J
    end

    # Solve
    # cache, meef, kepf, time, retcode = generate_qlaw_transfer(
    #     qLawPs, cost, QLawIndirectOptimization.ThreadedPSO; 
    #     max_time        = 2*3600.0, 
    #     show_trace      = true, 
    #     num_particles   = 30,
    # )

    # First solution found (50 particles with only 3600 seconds run time)
    qLawPs.oeW .= [6.714806348309516,7.957793270469378,5.25926959747461,1.9088836107960239,0.0]
    qLawPs.ηr = 0.09134109317011999
    cache, meef, kepf, time, retcode = generate_qlaw_transfer(qLawPs)

    # Save solution information
    jldsave(
        joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLike_MinFuel.jld2");
        cache = cache, params = qLawPs,
    )
    dump_to_mat(cache, joinpath(@__DIR__, "..", "data", "TAES", "mat", "MolniyaLike_MinFuel.mat"))

    # Generate figures
    plot_transfer(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "MolniyaLike_MinFuel_xy.png"), cache, qLawPs; 
        axes = [1,2], linewidth=0.2,
    )
    plot_transfer(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "MolniyaLike_MinFuel_xz.png"), cache, qLawPs; 
        axes = [1,3], linewidth=0.2,
    )
    plot_transfer(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "MolniyaLike_MinFuel_yz.png"), cache, qLawPs; 
        axes = [2,3], linewidth=0.2,
    )

    @infiltrate
    return (cache, qLawPs)
end

cache, ps = main()