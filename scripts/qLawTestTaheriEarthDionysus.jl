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
oeW          = [1.0, 1.0, 1.0, 1.0, 1.0] 

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
                ηr_tol      = 0.4,
                meeParams   = meeParams,
                spaceCraft  = spaceCraft,
                desolver    = Vern7(),
                maxRevs     = 400.0,
                integStep   = 5.0,
                writeData   = true,
                type        = :SD)

# # Define upper and lower bounds
# LB          = [0.1,   0.1,  0.1,  0.1,  0.1, 0.0]
# UB          = [10.0, 10.0, 10.0, 10.0, 10.0, 0.9]

# # Define cost function 
# function psoCost(x, ps_in)
#     # Create copy to remain thread safe
#     ps = deepcopy(ps_in)

#     # Set variables
#     ps.oeW[1:5] .= @view x[1:5]
#     ps.ηr = x[6]

#     # Run sim
#     J = try
#         tf, kepf, retcode = qLaw(qLawPs)

#         # Compute errors
#         aerr = ps.Ws[1]*abs(kepf[1] / ps.meePs.LU - ps.oet[1])
#         eerr = ps.Ws[2]*abs(kepf[2] - ps.oet[2])
#         ierr = ps.Ws[3]*abs(kepf[3] * pi/180.0 - ps.oet[3])
#         Ωerr = ps.Ws[4]*abs(acos(cos(kepf[4] * pi/180.0 - ps.oet[4])))
#         ωerr = ps.Ws[5]*abs(acos(cos(kepf[5] * pi/180.0 - ps.oet[5])))
#         errSum = aerr + eerr + ierr + Ωerr + ωerr

#         # Minimize time
#         J   = tf

#         # Penalize for failing
#         if retcode != :success
#             J += 1000.0 * errSum
#         end
#         J
#     catch
#         1e9
#     end

#     return J
# end

# # Perform PSO optimization
# prob = Problem(x -> psoCost(x,qLawPs), LB, UB)
# opts = Options(display = true, maxIters = 1000, maxStallIters = 50, 
#             funcTol = 1e-6, useParallel = true)
# pso  = PSO(prob; numParticles = 200)
# res  = optimize!(pso, opts)

# # Run simulation with final solution
# qLawPs.oeW .= res.xbest[1:5]
# qLawPs.ηr   = res.xbest[6]
# qLawPs.writeDataToFile = true

# Run QLaw sim
tf, kepf, retcode = qLawOriginal(qLawPs)

