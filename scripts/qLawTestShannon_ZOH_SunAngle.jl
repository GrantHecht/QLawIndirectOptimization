using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
using Heuristics
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2000-03-22T00:00:00")

# Compute thrust used by Shannon et al.
P    = 5.0*1000.0   # [W]
Isp  = 1800.0       # [s]
g0   = 9.80664      # [m/s^2]
η    = 0.55 
tMax = 2*η*P / (g0 * Isp)

# Define parameters for EOMs
μs          = 3.986e5
m0          = 1200.0
ephemDays   = 1000.0
ephemTspan  = (initEpoch - 100.0, initEpoch + ephemDays*86400.0)
nPoints     = ceil(Int64, 2*ephemDays)
tbEphems    = Ephemerides(ephemTspan, nPoints, [10,301], 399, "J2000")
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs,
                        thirdBodyEphemerides = tbEphems)
spaceCraft  = SimpleSpacecraft(m0, m0, tMax, Isp)

# Define initial and target orbital elements
mee0        = SVector(11359.07, 0.7306, 0.0, 0.2539676, 0.0, 0.0)
kep0, f     = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Keplerian, μs)
kept        = [42165.0, 0.01, 0.01, 0.0, 0.0]

# Convert angles in initial kep state to deg
kep0d        = Vector(kep0)
kep0d[3:6] .*= 180.0 / pi

# Define error weights
oeW         = [1.193, 2.402, 8.999, 0.0, 0.0] 
#oeW         = [10.0, 
#               4.092303319221676, 
#               4.723336735448565, 
#               4.546948268457167, 
#               2.8497164860046915]
#ηr          = 0.5922481336132026

# Define tolerance on targeted elements
atol        = 20.0
etol        = 0.001
itol        = 0.01
Ωtol        = 0.01
ωtol        = 0.01
tolVec      = [atol,etol,itol,Ωtol,ωtol]

# Construct qLaw parameters
qLawPs       = qLawParams(kep0d, kept;
                oeW         = oeW,
                oeTols      = tolVec,
                ηr_tol      = 0.0151,
                meeParams   = meeParams,
                spaceCraft  = spaceCraft,
                desolver    = Vern7(),
                maxRevs     = 500.0,
                integStep   = 5.0,
                writeData   = true,
                type        = :QDSAA,
                eSteps      = 10,
                eclipsing   = true,
                thrustSunAngleConstraint = true,
                thrustSunAngle = 60.0*pi/180.0,
                onlyWriteDataAtSteps = true)

# Run QLaw sim
#tf, kepf, retcode = qLawOriginal(qLawPs)
tf, kepf, retcode = qLawOriginal(qLawPs, :minfuel)
