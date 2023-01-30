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
tMax        = 1.0
Isp         = 3100.0
m0          = 300.0
meeParams   = MEEParams(initEpoch; LU = 384400.0, MU = 1.0, TU = 24.0*3600.0, μ = μs)
spaceCraft  = SimpleSpacecraft(m0, m0, tMax, Isp)

# Define initial and target orbital elements
μ           = AstroEOMs.getScaledGravityParameter(meeParams)
kep0        = SVector(7000.0 / meeParams.LU, 0.01, 0.05*pi/180, 0.0, 0.0, 0.0)
mee0        = AstroUtils.convertState(kep0, AstroUtils.Keplerian, AstroUtils.MEE, μ)
cart0       = AstroUtils.convertState(mee0, AstroUtils.MEE, AstroUtils.Cartesian, μ)
fullState0  = SVector(mee0[1], mee0[2], mee0[3], mee0[4], mee0[5], mee0[6], spaceCraft.initMass)
kept        = [42000.0 / meeParams.LU, 0.01, 10*pi/180, 0.0, 0.0]

# Define qLaw parameters
oeW          = [1.0, 1.0, 0.0, 0.0, 0.0] 
qLawPs       = qLawParams(kept, oeW, 0.0, 6578.0 / meeParams.LU, 1.0, μ,
                spaceCraft.tMax * meeParams.TU^2 / (1000.0*meeParams.MU*meeParams.LU),
                0.861, 360)

# Define EOMs
function qLawEOMs(u, p, t)
    # Grab parameters
    meeParams   = p[1]
    spaceCraft  = p[2]
    qLawPs      = p[3]

    # Construct MEE state
    mee         = SVector(u[1],u[2],u[3],u[4],u[5],t)

    # Compute acceleration due to thrust
    if qLawPs.coasting == true
        at   = SVector(0.0,0.0,0.0)
        umag = 0.0
    else
        α  = qLawPs.α
        β  = qLawPs.β
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
scb  = SavingCallback(save_func,saved_values)

# Deal with Sundman transformation
fullState0s = SVector(fullState0[1], fullState0[2], fullState0[3], 
                fullState0[4], fullState0[5], 0.0, fullState0[7])

# Specify independant variable span and info
nRevs       = 600.0
step        = 5.0 * pi / 180
Lspan       = (fullState0[6], fullState0[6] + nRevs*2*pi)

# Allocate storage arrays
n           = ceil(Int,90.0*nRevs)
Ls          = range(Lspan[1], Lspan[2]; length = n)
cart        = fill(NaN, n, 7)
mee         = fill(NaN, n, 7)
kep         = fill(NaN, n, 7)
coast       = fill(true, n)

# Begin integration loop
done = false
idx  = 1
L0   = Lspan[1]
Lf   = L0 + step
while !done
    # Compute kep state
    mee     = SVector(fullState0[1], fullState0[2], fullState0[3], 
                    fullState0[4], fullState0[5], L0)
    cart    = AstroUtils.convertState(mee, AstroUtils.Mee,
                    AstroUtils.Cartesian, qLawPs.μ)
    kep, fl = AstroUtils.convertState(cart, AstroUtils.Cartesian,
                    AstroUtils.Keplerian, qLawPs.μ)

    # Compute coasting
    qLawPs.coasting = false
    if qLawPs.ηr > 0.0
        qLawPs.coasting = 
            qLawCoastContinuousCallbackCheck(kep, fullState0[7], qLawPs) > 0.0 ? false : true
    end

    # Compute thrust angles
    if !coasting
        α, β = qLawThrustAngles(kep[1], kep[2], kep[3], kep[4], kep[5], kep[6], 
                    fullState0[7], qLawPs)
        qLawPs.α = α
        qLawPs.β = β
    end

    # Perform numerical integration
    prob = ODEProblem(qLawEOMs, fullState0s, Lspan, (meeParams,spaceCraft,qLawPs))
    sol  = solve(prob, Vern9(), reltol=1e-12, abstol=1e-12, 
        dtmax = 20*pi/180, callback = scb)

    # Save info
    while Ls[idx] <= sol.t[end]
        mees_us     = sol(Ls[idx])
        mee_us      = SVector(mees_us[1], mees_us[2], mees_us[3],
                        mees_us[4], mees_us[5], Ls[idx], mees_us[7])

        mee[idx,1]      = meeParams.LU * mee_us[1]
        mee[idx,2:7]   .= mee_us[2:7] 

        cart_us         = AstroUtils.convertState(mee_us, AstroUtils.MEE, AstroUtils.Cartesian, μ)
        cart[idx,1:3]  .= meeParams.LU*view(cart_us,1:3)
        cart[idx,4:6]  .= meeParams.LU*view(cart_us,4:6)/meeParams.TU
        cart[idx,7]     = mee_us[7]

        kep_us,fff      = AstroUtils.convertState(cart_us, AstroUtils.Cartesian, AstroUtils.Keplerian, μ)
        kep[idx,1]      = meeParams.LU*kep_us[1]
        kep[idx,2:6]   .= view(kep_us, 2:6)
        kep[idx,7]      = mee_us[7]

        # Deal with coasting
        idxFound = false
        j        = 0
        while !idxFound 
            j += 1
            if Ls[idx] >= saved_values.t[j] && Ls[idx] <= saved_values.t[j + 1]
                idxFound = true
            end
        end
        coast[idx]      = saved_values.saveval[j]

        # Increment index
        idx += 1
    end

    # Check if we need to stop
    if Lf >= Lspan[2]
        done = true
end

# Write data to files
open(datadir("kep.txt"),   "w") do io; writedlm(io,   kep); end
open(datadir("mee.txt"),   "w") do io; writedlm(io,   mee); end
open(datadir("cart.txt"),  "w") do io; writedlm(io,  cart); end
open(datadir("coast.txt"), "w") do io; writedlm(io, coast); end

#plot(cart[:,1],cart[:,2],cart[:,3])
#plot(ts,kep[:,1])
