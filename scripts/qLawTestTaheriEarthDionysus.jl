using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
using Heuristics
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2005-01-30T00:00:00")

# Define parameters for EOMs
μs          = 132712440018.0
tMax        = 0.32
Isp         = 3000.0
m0          = 4000.0
mp          = 3000.0
LU          = 1.496e8
TU          = 365.0*24.0*3600.0 / (2*pi)
meeParams   = MEEParams(initEpoch; LU = LU, MU = 1.0, TU = TU, μ = μs)
spaceCraft  = SimpleSpacecraft(m0, mp, tMax, Isp)

# Define initial and target orbital elements
cart0       = SVector(-3637871.081, 147099798.784, -2261.441,
                -30.265097, -0.8486854, 0.0000505)
kep0,f      = AstroUtils.convertState(cart0, AstroUtils.Cartesian, AstroUtils.Keplerian, μs)
kep0d       = Vector(kep0)
kep0d[3:6].*= 180.0 / pi
kept        = [2.2*LU, 0.542, 13.6, 82.2, 204.2]

# Define qLaw parameters
oeW          = [1.0, 1.0, 1.0, 1.0, 50.0] 

# Define tolerance on targeted elements
atol        = 1000.0
etol        = 0.001
itol        = 0.01
Ωtol        = 0.01
ωtol        = 0.01
tolVec      = [atol,etol,itol,Ωtol,ωtol]

# Construct qLaw parameters
qLawPs       = qLawParams(copy(kep0d), copy(kept);
                oeW         = oeW,
                oeTols      = tolVec,
                ηr_tol      = 0.2,
                ηa_tol      = 0.0,
                meeParams   = meeParams,
                spaceCraft  = spaceCraft,
                desolver    = Vern7(),
                maxRevs     = 100.0,
                integStep   = 5.0,
                writeData   = true,
                type        = :SD)

# Run QLaw sim
tf, kepf, retcode = qLawOriginal(qLawPs)

