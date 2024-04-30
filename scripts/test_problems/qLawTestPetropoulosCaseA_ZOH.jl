using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
furnshDefaults()

function main()
    # Compute initial epoch
    initEpoch   = utc2et("2022-10-07T12:00:00")

    # Define parameters for EOMs
    μs          = 3.986e5
    tMax        = 1.0
    Isp         = 3100.0
    m0          = 300.0
    mp          = 290.0
    meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs)
    spaceCraft  = SimpleSpacecraft(m0, mp, tMax, Isp)

    # Define initial and target orbital elements
    kep0        = [7000.0,  0.01, 0.05, 0.0, 0.0, 0.0]
    kept        = [42000.0, 0.01, 0.0, 0.0, 0.0]

    # Define qLaw parameters
    oeW         = [1.0, 1.0, 0.0, 0.0, 0.0]

    # Define tolerance on targeted elements
    atol        = 10.0
    etol        = 0.01
    itol        = 0.01
    Ωtol        = 0.01
    ωtol        = 0.01
    tolVec      = [atol,etol,itol,Ωtol,ωtol]

    # Define qLaw parameters
    qLawPs       = qLawParams(kep0, kept;
                    oeW         = oeW, 
                    oeTols      = tolVec,
                    ηr_tol      = 0.167,
                    ηa_tol      = 0.0,
                    meeParams   = meeParams,
                    spaceCraft  = spaceCraft,
                    desolver    = Vern7(),
                    maxRevs     = 200.0,
                    integStep   = 0.1,
                    writeData   = true,
                    type        = :QDUC)

    # Run QLaw sim 
    tf, kepf, retcode = qLawOriginal(qLawPs)
end

main()