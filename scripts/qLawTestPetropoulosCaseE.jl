using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays
using DifferentialEquations, DiffEqCallbacks, Plots
using DelimitedFiles
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2022-10-07T12:00:00")

# Include source
include(srcdir("includeSource.jl"))

# Define parameters for EOMs
μs          = 3.986e5
tMax        = 2.0
Isp         = 2000.0
m0          = 2000.0
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs)
spaceCraft  = SimpleSpacecraft(m0, m0, tMax, Isp)

# Define initial and target orbital elements
μ           = AstroEOMs.getScaledGravityParameter(meeParams)
kep0        = SVector(24505.9 / meeParams.LU, 0.725, 0.06*pi/180, 0.0, 0.0, 0.0)
mee0        = AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μ)
cart0       = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Cartesian, μ)
fullState0  = SVector(mee0[1], mee0[2], mee0[3], mee0[4], mee0[5], mee0[6], spaceCraft.initMass)
kept        = [26500.0 / meeParams.LU, 0.7, 116*pi/180, 270.0*pi/180, pi]

# Define qLaw parameters
oeW          = [1.0, 1.0, 2.0, 1.0, 1.0] 
qLawPs       = qLawParams(kept, oeW, 1.0, 6578.0 / meeParams.LU, 1.0, μ,
                spaceCraft.tMax * meeParams.TU^2 / (1000.0*meeParams.MU*meeParams.LU),
                0.652, 360)

# Define EOMs
function qLawEOMs(u, p, t)
    # Grab parameters
    meeParams   = p[1]
    spaceCraft  = p[2]
    qLawPs      = p[3]

    # Construct MEE state
    mee         = SVector(u[1],u[2],u[3],u[4],u[5],t)

    # Compute keplerian elements
    cart        = AstroUtils.convertState(mee, AstroUtils.MEE, 
                        AstroUtils.Cartesian, qLawPs.μ)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, qLawPs.μ)

    # Now basing coasting control on callback
    if qLawPs.coasting == true
        at   = SVector(0.0,0.0,0.0)
        umag = 0.0
    else
        α,β  = qLawThrustAngles(kep[1],kep[2],kep[3],kep[4],kep[5],kep[6],u[7],qLawPs)
        at = SVector(qLawPs.tMax*cos(β)*sin(α) / u[7],
                     qLawPs.tMax*cos(β)*cos(α) / u[7],
                     qLawPs.tMax*sin(β) / u[7])
        umag = qLawPs.tMax 
    end

    # Compute state dynamics
    dmee    = AstroEOMs.meeEomControl(mee,meeParams,u[6],at)
    dLinv   = 1.0 / dmee[6]

    # Compute sundman transformed dynamics
    dmees   = SVector(dmee[1]*dLinv, dmee[2]*dLinv, dmee[3]*dLinv,
                dmee[4]*dLinv, dmee[5]*dLinv, dLinv)

    # Compute mass dynamics
    cs  = spaceCraft.c * meeParams.TU / (1000.0 * meeParams.LU)
    dm  = -umag / cs
    dms = dm*dLinv

    # Return full state dynamics
    return SVector(dmees[1],dmees[2],dmees[3],dmees[4],dmees[5],dmees[6],dms)
end

# Save function
saved_values = SavedValues(Float64,Bool)
function save_func(u,t,integrator)
    return (integrator.p[3].coasting)
end

# Define integration termination callback
function term_condition(u,t,integrator)
    ps  = integrator.p
    meeParams   = ps[1]
    spaceCraft  = ps[2]
    qLawPs      = ps[3]

    # Construct mee state
    mee         = SVector(u[1],u[2],u[3],u[4],u[5],t)

    # Compute keplerian elements
    cart        = AstroUtils.convertState(mee, AstroUtils.MEE, 
                        AstroUtils.Cartesian, qLawPs.μ)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, qLawPs.μ)

    # Set tolerances
    atol        = 10.0 / meeParams.LU
    etol        = 0.001
    itol        = 0.01*pi / 180
    Ωtol        = 0.01*pi / 180
    ωtol        = 0.01*pi / 180
    tolVec      = SVector(atol,etol,itol,Ωtol,ωtol)

    # Set weights
    Wa          = qLawPs.oeW[1] > 0.0 ? 1.0 : 0.0
    We          = qLawPs.oeW[2] > 0.0 ? 1.0 : 0.0
    Wi          = qLawPs.oeW[3] > 0.0 ? 1.0 : 0.0
    WΩ          = qLawPs.oeW[4] > 0.0 ? 1.0 : 0.0
    Wω          = qLawPs.oeW[5] > 0.0 ? 1.0 : 0.0
    Ws          = SVector(Wa,We,Wi,WΩ,Wω)

    # Compute return value
    val         = maximum(Ws .* (abs.(kep[1:5] - qLawPs.oet) - tolVec))
    print(string(meeParams.LU*(abs(kep[1] - qLawPs.oet[1]))) * "\t" *
          string(abs(kep[2] - qLawPs.oet[2])) * "\t" * 
          string(abs(kep[3] - qLawPs.oet[3])) * "\t" * 
          string(t / (2*pi)) * "\n")
    #print(string(t) * "\n")
    return val
end

term_affect!(integrator) = terminate!(integrator)

# Define qlaw coast switching callback
function coast_condition(u,t,integrator)
    ps  = integrator.p
    meeParams   = ps[1]
    spaceCraft  = ps[2]
    qLawPs      = ps[3]

    # Construct mee state
    mee     = SVector(u[1],u[2],u[3],u[4],u[5],t)

    if qLawPs.ηr > 0.0
        # Compute keplerian elements
        cart        = AstroUtils.convertState(mee, AstroUtils.MEE, 
                            AstroUtils.Cartesian, qLawPs.μ)
        kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                            AstroUtils.Keplerian, qLawPs.μ)

        # Check effectivity
        val = qLawCoastContinuousCallbackCheck(kep,u[7],qLawPs)
    else
        return 1.0
    end
end

function coast_affect_pos!(integrator)
    integrator.p[3].coasting = false
    return nothing
end

function coast_affect_neg!(integrator)
    integrator.p[3].coasting = true
    return nothing
end

# Construct Callbacks
tcb  = ContinuousCallback(term_condition,term_affect!)
ccb  = ContinuousCallback(coast_condition,coast_affect_pos!,coast_affect_neg!)
scb  = SavingCallback(save_func,saved_values)

# Set initial coasting flag
val  = qLawCoastContinuousCallbackCheck(kep0,fullState0[7],qLawPs)
if qLawPs.ηr > 0.0
    if val > 0 
        qLawPs.coasting = false
    else
        qLawPs.coasting = true
    end
else
    qLawPs.coasting = false
end

# Deal with Sundman transformation
fullState0s = SVector(fullState0[1], fullState0[2], fullState0[3], 
                fullState0[4], fullState0[5], 0.0, fullState0[7])
nRevs       = 300.0
Lspan       = (fullState0[6], fullState0[6] + nRevs*2*pi)

# Perform numerical integration
prob = ODEProblem(qLawEOMs, fullState0s, Lspan, (meeParams,spaceCraft,qLawPs))
sol  = solve(prob, Vern9(), reltol=1e-8, abstol=1e-8, 
    dtmax = 20*pi/180, callback = CallbackSet(tcb,ccb,scb))

# Grab relavent info
ts   = range(0.0, sol.t[end]; length = floor(Int, 90.0*sol.t[end]/(2*pi)))
cart = zeros(length(ts),7)
mee  = zeros(length(ts),7)
kep  = zeros(length(ts),7)
coast= zeros(length(ts))
for i in eachindex(ts)
    mees_us     = sol(ts[i])
    mee_us      = SVector(mees_us[1], mees_us[2], mees_us[3],
                    mees_us[4], mees_us[5], ts[i], mees_us[7])

    mee[i,1]    = meeParams.LU * mee_us[1]
    mee[i,2:7]  .= mee_us[2:7] 

    cart_us     = AstroUtils.convertState(mee_us, AstroUtils.MEE, AstroUtils.Cartesian, μ)
    cart[i,1:3] .= meeParams.LU*view(cart_us,1:3)
    cart[i,4:6] .= meeParams.LU*view(cart_us,4:6)/meeParams.TU
    cart[i,7]   = mee_us[7]

    kep_us,fff  = AstroUtils.convertState(cart_us, AstroUtils.Cartesian, AstroUtils.Keplerian, μ)
    kep[i,1]    = meeParams.LU*kep_us[1]
    kep[i,2:6]  .= view(kep_us, 2:6)
    kep[i,7]    = mee_us[7]

    # Deal with coasting
    idxFound = false
    j        = 0
    while !idxFound 
        j += 1
        if ts[i] >= saved_values.t[j] && ts[i] <= saved_values.t[j + 1]
            idxFound = true
        end
    end
    coast[i] = saved_values.saveval[j][1]
end
#coast = [saved_values.saveval[i][1] for i in eachindex(saved_values.saveval)]

# Write data to files
open(datadir("kep.txt"),   "w") do io; writedlm(io,   kep); end
open(datadir("mee.txt"),   "w") do io; writedlm(io,   mee); end
open(datadir("cart.txt"),  "w") do io; writedlm(io,  cart); end
open(datadir("coast.txt"), "w") do io; writedlm(io, coast); end

#plot(cart[:,1],cart[:,2],cart[:,3])
#plot(ts,kep[:,1])
