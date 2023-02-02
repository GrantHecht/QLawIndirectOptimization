using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
using QLawIndirectOptimization
using Heuristics
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2022-10-07T12:00:00")

# Compute thrust used by Shannon et al.
P    = 5.0*1000.0   # [W]
Isp  = 1800.0       # [s]
g0   = 9.80664      # [m/s^2]
η    = 0.55 
tMax = 2*η*P / (g0 * Isp)

# Define parameters for EOMs
μs          = 3.986e5
m0          = 1200.0
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs)
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
                ηr_tol      = 0.01,
                meeParams   = meeParams,
                spaceCraft  = spaceCraft,
                desolver    = Vern7(),
                maxRevs     = 1000.0,
                integStep   = 5.0,
                writeData   = true)

# Define upper and lower bounds
LB          = [0.1,   0.1,  0.1,  0.1,  0.1, 0.0]
UB          = [10.0, 10.0, 10.0, 10.0, 10.0, 0.9]

# Define cost function 
function psoCost(x, ps_in)
    # Create copy to remain thread safe
    ps = deepcopy(ps_in)

    # Set variables
    ps.oeW[1:5] .= @view x[1:5]
    ps.ηr = x[6]

    # Run sim
    J = try
        tf, kepf, retcode = qLaw(qLawPs)

        # Compute errors
        aerr = ps.Ws[1]*abs(kepf[1] / ps.meePs.LU - ps.oet[1])
        eerr = ps.Ws[2]*abs(kepf[2] - ps.oet[2])
        ierr = ps.Ws[3]*abs(kepf[3] * pi/180.0 - ps.oet[3])
        Ωerr = ps.Ws[4]*abs(acos(cos(kepf[4] * pi/180.0 - ps.oet[4])))
        ωerr = ps.Ws[5]*abs(acos(cos(kepf[5] * pi/180.0 - ps.oet[5])))
        errSum = aerr + eerr + ierr + 1000.0*Ωerr + 1000.0*ωerr

        # Minimize time
        J   = -kepf[7]

        # Penalize for failing
        if retcode != :success
            J += 1000.0 * errSum
        end
        J
    catch
        1e9
    end

    return J
end

# Perform PSO optimization
prob = Problem(x -> psoCost(x,qLawPs), LB, UB)
opts = Options(display = true, maxIters = 1000, maxStallIters = 50, 
            funcTol = 1e-6, useParallel = true)
pso  = PSO(prob; numParticles = 64)
res  = optimize!(pso, opts)

# Run simulation with final solution
qLawPs.oeW .= res.xbest[1:5]
qLawPs.ηr   = res.xbest[6]
qLawPs.writeDataToFile = true

# Run QLaw sim
tf, kepf, retcode = qLaw(qLawPs)
