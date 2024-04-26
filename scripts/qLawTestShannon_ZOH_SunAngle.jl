using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations
using QLawIndirectOptimization
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
meeParams   = MEEParams(
    initEpoch; 
    LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, 
    μ = μs,
    thirdBodyEphemerides = tbEphems, 
    nonsphericalGravity = true,
)
spaceCraft  = SimpleSpacecraft(m0, m0, tMax, Isp)

# Define initial and target orbital elements
mee0        = SVector(11359.07, 0.7306, 0.0, 0.2539676, 0.0, 0.0)
kep0        = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Keplerian, μs)
kept        = [42165.0, 0.01, 0.01, 0.0, 0.0]

# Convert angles in initial kep state to deg
kep0d        = Vector{Float64}(kep0)
kep0d[3:6] .*= 180.0 / pi

# Define error weights
#oeW         = [1.193, 2.402, 8.999, 0.0, 0.0] 

# Minimum fuel PSO solution
#oeW         = [10.0,
#               6.159,
#               9.65, 
#               0.0, 0.0]
#ηr          = 0.374423

# Minimum time PSO solution (old)
#oeW         = [9.889426673974159,
#               4.284182682658705,
#               4.7885410198537945,
#               0.0,
#               0.0]
#ηr          = 0.0

# Minimum time PSO solution
oeW         = [3.6648209940624583,
               1.601742021805308,
               9.893525941369818, 0.0, 0.0]
ηr          = 0.2

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
    oeW                      = oeW,
    oeTols                   = tolVec,
    ηr_tol                   = 0.1,
    meeParams                = meeParams,
    spaceCraft               = spaceCraft,
    desolver                 = Vern7(),
    maxRevs                  = 800.0,
    integStepOpt             = 10.0,
    integStepGen             = 0.1,
    writeData                = true,
    type                     = :QDSAA,
    eSteps                   = 10,
    eclipsing                = true,
    thrustSunAngleConstraint = true,
    thrustSunAngle           = 50.0*pi/180.0,
    onlyWriteDataAtSteps     = true,
    savedStatesAtSteps       = 10,
)

# Run QLaw sim
#meefs, kepfus, tf, retcode = qLawOriginal(qLawPs)

#meefs_fast, kepfus_fast, tf_fast, retcode_fast = generate_qlaw_transfer(qLawPs)

#cache = QLawIndirectOptimization.QLawTransferCache()
#meefs_fast, kepfus_fast, tf_fast, retcode_fast = generate_qlaw_transfer(qLawPs, cache)

#tf, kepf, retcode = qLawOriginal(qLawPs, :mintime; maxTime = 5.0*3600.0)

# Define cost
function cost(state, time, retcode)
    J = time - 5*state[7]
    if retcode != :success
        J += 5000.0
    end
    return J
end

# Solve
cache, meef, kepf, time, retcode = generate_qlaw_transfer(
    qLawPs, cost; 
    max_time        = 60.0, 
    show_trace      = true, 
    num_particles   = 50,
)
fig = plot_transfer(cache, qLawPs)