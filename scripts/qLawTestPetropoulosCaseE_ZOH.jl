using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2022-10-07T12:00:00")

# Define parameters for EOMs
μs          = 3.986e5
tMax        = 2.0
Isp         = 3100.0
m0          = 3000.0
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs)
spaceCraft  = SimpleSpacecraft(m0, 0.8*m0, tMax, Isp)

# Define initial and target orbital elements
μ           = AstroEOMs.getScaledGravityParameter(meeParams)
kep0        = [24505.9, 0.725, 0.06, 0.0, 0.0, 0.0]
mee0        = AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μs)
kept        = [26500.0, 0.7, 116.0, 175.0, 270.0]

# Define qLaw parameters
oeW          = [1.0, 1.0, 1.0, 1.0, 100.0] 

# Define tolerance on targeted elements
atol        = 1000.0
etol        = 0.001
itol        = 0.01
Ωtol        = 0.01
ωtol        = 0.01
tolVec      = [atol,etol,itol,Ωtol,ωtol]

# Construct qLaw parameters
qLawPs       = qLawParams(copy(kep0), copy(kept);
                oeW         = oeW,
                Wp          = 1.0,
                oeTols      = tolVec,
                ηr_tol      = 0.1,
                ηa_tol      = 0.0,
                meeParams   = meeParams,
                spaceCraft  = spaceCraft,
                desolver    = Vern7(),
                maxRevs     = 5000.0,
                integStep   = 1.0,
                writeData   = true,
                type        = :SD)

tf, kepf, retcode = qLawOriginal(qLawPs)