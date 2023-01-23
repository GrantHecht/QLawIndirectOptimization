using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays, DifferentialEquations
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2022-10-07T12:00:00")

# Include source
include(srcdir("includeSource.jl"))

# Define parameters for EOMs
μ           = 3.986e5
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μ)
spaceCraft  = SimpleSpacecraft(300.0, 300.0, 1.0, 3100.0)

# Define initial and target orbital elements
kep0        = SVector(7000.0,0.01,0.05*pi/180,0.0,0.0,0.0)
mee0        = AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μ)
fullState0  = SVector(mee0[1],mee0[2],mee0[3],mee0[4],mee0[5],mee0[6],spaceCraft.initMass)
kept        = [42000.0,0.01,0.05*pi/180,0.0,0.0]

# Define qLaw parameters
oeW         = [1.0,1.0,1.,1.,1.0] 
qLawPs      = qLawParams(kept,oeW,AstroEOMs.getScaledGravityParameter(meeParams),
                spaceCraft.tMax * meeParams.TU / (1000.0*meeParams.MU*meeParams.LU))
qLawPs.rpmin = 6300.0

# Define EOMs
function qLawEOMs(u, p, t)
    # Grab parameters
    meeParams   = p[1]
    spaceCraft  = p[2]
    qLawPs      = p[3]

    # Compute keplerian elements
    cart        = AstroUtils.convertState(u, AstroUtils.MEE, 
                        AstroUtils.Cartesian, meeParams.mu)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, meeParams.mu)

    # Compute qLaw control
    α,β,coast = qLaw(kep,qLawPs)

    # Compute acceleration vector
    if coast
        au   = SVector(0.0,0.0,0.0)
        umag = 0.0
    else
        au = SVector(qLawPs.tMax*cos(β)*sin(α) / u[7],
                     qLawPs.tMax*cos(β)*cos(α) / u[7],
                     qLawPs.tMax*sin(β) / u[7])
        umag = qLawPs.tMax 
    end

    # Compute state dynamics
    dmee = AstroEOMs.meeEomControl(u,meeParams,t,au)

    # Compute mass dynamics
    dm  = -umag*meeParams.TU / (spaceCraft.c * meeParams.MU) 

    # Return full state dynamics
    return SVector(dmee[1],dmee[2],dmee[3],dmee[4],dmee[5],dmee[6],dm)
end

# Perform numerical integration
prob = ODEProblem(qLawEOMs, fullState0, (0.0, 25.0), (meeParams,spaceCraft,qLawPs))
sol  = solve(prob)