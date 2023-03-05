using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2012-12-23T12:00:00")

# Define parameters for EOMs
μs          = 3.986e5
tMax        = 0.5
Isp         = 3100.0
m0          = 100.0
mp          = 90.0
LU          = 6378.0
TU          = sqrt(LU^3 / μs)
meeParams   = MEEParams(initEpoch; LU = LU, MU = 1.0, TU = TU, μ = μs, 
                costFunction = AstroEOMs.MinimumFuelMayer)
spaceCraft  = SimpleSpacecraft(m0, mp, tMax, Isp)
sw          = SwitchingStruct(1.0, 1)

# Define initial and target orbital elements
cart0       = SVector(6378.9, 0.0, 0.0,
                0.0, 10.0258, 1.231)
kep0,f      = AstroUtils.convertState(cart0, AstroUtils.Cartesian, AstroUtils.Keplerian, μs)
mee0        = Vector(AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μs))
kep0d       = Vector(kep0)
for i in 3:6
    kep0d[i] = mod(kep0d[i], 2.0*pi) * 180.0 / pi
end
kept        = [42165.0, 0.01, 0.01, 0.0, 0.0]

# Define qLaw parameters
oeW          = [1.0, 1.0, 5.0, 0.0, 0.0] 

# Define tolerance on targeted elements
atol        = 20.0
etol        = 0.001
itol        = 0.01
Ωtol        = 0.01
ωtol        = 0.01
tolVec      = [atol,etol,itol,Ωtol,ωtol]

# Construct qLaw parameters
qLawPs       = qLawParams(copy(kep0d), copy(kept);
                oeW         = oeW,
                oeTols      = tolVec,
                ηr_tol      = 0.15,
                ηa_tol      = 0.0,
                meeParams   = meeParams,
                spaceCraft  = spaceCraft,
                desolver    = Vern7(),
                maxRevs     = 20.0,
                integStep   = 2.0*pi / 180.0,
                writeData   = true,
                returnData  = true,
                type        = :QDUC)

# Run QLaw sim
meeAtSteps, kepf, retcode = qLawOriginal(qLawPs)

# Scale initial mee state
mee0[1] /= meeParams.LU

# Construct interplants starting with state at ΔL = 2*pi
L0      = meeAtSteps[1,6]
istart  = 0
for i in axes(meeAtSteps,1)
    if abs(meeAtSteps[i,6] - L0 - 12*pi) < 1e-3
        global istart = i
        break
    end
end
interp  = FinalStateInterpolant(meeAtSteps, istart)

# Solve short problem with κ = 0.0
target      = getFinalStateAndTime(interp, 0.0)
λ0, retcode = minFuelMayerSolve(mee0, target[1:5], (0.0, target[6]), (spaceCraft,meeParams,sw))

# Perform PC continuation
pcContinuation(λ0, 1e-1, mee0, interp, (spaceCraft,meeParams,sw))