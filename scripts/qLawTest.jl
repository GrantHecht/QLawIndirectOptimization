using DrWatson
@quickactivate "QLawIndirectOptimization"

using AstroEOMs, AstroUtils, SPICE, StaticArrays, DifferentialEquations, DiffEqCallbacks, Plots
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2022-10-07T12:00:00")

# Include source
include(srcdir("includeSource.jl"))

# Define parameters for EOMs
μs          = 3.986e5
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs)
spaceCraft  = SimpleSpacecraft(1200.0, 1200.0, 0.3115799539654781, 1800.0)

# Define initial and target orbital elements
μ           = AstroEOMs.getScaledGravityParameter(meeParams)
#kep0        = SVector(24505.9 / meeParams.LU, 0.725, 0.06*pi/180, 0.0, 0.0, 0.0)
#mee0        = AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μ)
mee0        = SVector(11359.07 / meeParams.LU, 0.7306, 0.0, 0.2539676, 0.0, 0.0)
cart0       = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Cartesian, μ)
kep0, f     = AstroUtils.convertState(cart0, AstroUtils.Cartesian, AstroUtils.Keplerian, μ)
fullState0  = SVector(mee0[1], mee0[2], mee0[3], mee0[4], mee0[5], mee0[6], spaceCraft.initMass)
kept        = [42165.0 / meeParams.LU, 0.01, 0.01*pi/180, 270.0*pi/180, 180*pi/180]

# Define qLaw parameters
oeW          = [1.193, 2.402, 8.999, 0.0, 0.0] 
qLawPs       = qLawParams(kept, oeW, 1.0, 6578.0, 100.0, μ,
                spaceCraft.tMax * meeParams.TU^2 / (1000.0*meeParams.MU*meeParams.LU),
                0.0151, 360)

# Define EOMs
function qLawEOMs(u, p, t)
    # Grab parameters
    meeParams   = p[1]
    spaceCraft  = p[2]
    qLawPs      = p[3]

    # Compute keplerian elements
    cart        = AstroUtils.convertState(u, AstroUtils.MEE, 
                        AstroUtils.Cartesian, qLawPs.μ)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, qLawPs.μ)

    # Now basing coasting control on callback
    if qLawPs.coasting == true
        at   = SVector(0.0,0.0,0.0)
        umag = 0.0
    else
        au = qLawThrustUnitVector(kep[1],kep[2],kep[3],kep[4],kep[5],kep[6],qLawPs)
        at = SVector(qLawPs.tMax*au[1] / u[7],
                     qLawPs.tMax*au[2] / u[7],
                     qLawPs.tMax*au[3] / u[7])
        umag = qLawPs.tMax 
    end

    # Compute state dynamics
    dmee = AstroEOMs.meeEomControl(u,meeParams,t,at)

    # Compute mass dynamics
    cs  = spaceCraft.c * meeParams.TU / (1000.0 * meeParams.LU)
    dm  = -qLawPs.tMax / cs

    # Return full state dynamics
    return SVector(dmee[1],dmee[2],dmee[3],dmee[4],dmee[5],dmee[6],dm)
end

# Save function
saved_values = SavedValues(Float64,Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64})
function save_func(u,t,integrator)
    # Grab parameters
    meeParams   = integrator.p[1]
    spaceCraft  = integrator.p[2]
    qLawPs      = integrator.p[3]

    # Compute keplerian elements
    cart        = AstroUtils.convertState(u, AstroUtils.MEE, 
                        AstroUtils.Cartesian, qLawPs.μ)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, qLawPs.μ)

    # Compute qLaw control
    au,coast = qLaw(kep, qLawPs)

    return (coast, kep[1], kep[2], kep[3], kep[4], kep[5], kep[6])
end

# Define integration termination callback
function term_condition(u,t,integrator)
    ps  = integrator.p
    meeParams   = ps[1]
    spaceCraft  = ps[2]
    qLawPs      = ps[3]

    # Compute keplerian elements
    cart        = AstroUtils.convertState(u, AstroUtils.MEE, 
                        AstroUtils.Cartesian, qLawPs.μ)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, qLawPs.μ)

    # Set tolerances
    atol        = 10 / meeParams.LU
    etol        = 0.01
    itol        = 0.01*pi / 180
    Ωtol        = 0.01*pi / 180
    ωtol        = 0.01*pi / 180
    tolVec      = SVector(atol,etol,itol,Ωtol,ωtol)

    # Compute keplerian element errors (FOR DESCRETE CALLBACK)
    #err         = abs.(qLawPs.oeW.*(u[1:5] - qLawPs.oet))
    #err[1]     *= meeParams.LU

    #stop        = err[1] < 10 && err[2] < 0.01 && err[3] < 5*pi/180 && err[4] < 5*pi/180 && err[5] < 5*pi/180
    #return stop

    # Compute return value
    val         = sum(qLawPs.oeW .* (abs.(u[1:5] - qLawPs.oet) - tolVec))
    print(string(val) * "\n")
    return val
end

term_affect!(integrator) = terminate!(integrator)

# Define qlaw coast switching callback
function coast_condition(u,t,integrator)
    ps  = integrator.p
    meeParams   = ps[1]
    spaceCraft  = ps[2]
    qLawPs      = ps[3]

    # Compute keplerian elements
    cart        = AstroUtils.convertState(u, AstroUtils.MEE, 
                        AstroUtils.Cartesian, qLawPs.μ)
    kep,fl      = AstroUtils.convertState(cart, AstroUtils.Cartesian, 
                        AstroUtils.Keplerian, qLawPs.μ)

    # Check effectivity
    val = qLawCoastContinuousCallbackCheck(kep,qLawPs)
    return val
end

function coast_affect_pos!(integrator)
    integrator.p[3].coasting = false
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
val  = qLawCoastContinuousCallbackCheck(kep0,qLawPs)
if val > 0
    qLawPs.coasting = false
else
    qLawPs.coasting = true
end

# Perform numerical integration
prob = ODEProblem(qLawEOMs, fullState0, (0.0, 200.0), (meeParams,spaceCraft,qLawPs))
sol  = solve(prob, Vern9(), reltol=1e-12, abstol=1e-12, dtmax = 3600.0 / 86400.0,callback = CallbackSet(tcb,ccb,scb))

ts   = range(0.0, sol.t[end]; length = 10000)
cart = zeros(length(ts),6)
kep  = zeros(length(ts),6)
for i in eachindex(ts)
    cart_us     = AstroUtils.convertState(sol(ts[i]), AstroUtils.MEE, AstroUtils.Cartesian, μ)
    cart[i,1:3] .= meeParams.LU*view(cart_us,1:3)
    cart[i,4:6] .= meeParams.LU*view(cart_us,4:6)/meeParams.TU

    kep_us,fff  = AstroUtils.convertState(cart_us, AstroUtils.Cartesian, AstroUtils.Keplerian, μ)
    kep[i,1]    = meeParams.LU*kep_us[1]
    kep[i,2:6]  .= view(kep_us, 2:6)
end

#plot(cart[:,1],cart[:,2],cart[:,3])
plot(ts,kep[:,1])