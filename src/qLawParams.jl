# Define some qlaw parameters
mutable struct qLawParams
    # Initial orbital elements
    oei::Vector{Float64}

    # Target orbital elements
    oet::Vector{Float64}

    # Target orbital element switch
    oets::Vector{Bool}

    # Convergence tolerances
    oetols::Vector{Float64}

    # QLaw parameters
    oeW::Vector{Float64}
    Wp::Float64
    rpmin::Float64
    k::Float64
    m_petro::Float64
    n_petro::Float64
    r_petro::Float64
    b_petro::Float64

    # Effectivity settings
    effectivity::Bool
    abs_effectivity::Bool
    rel_effectivity::Bool
    η_rel_tol::Bool
    η_abs_tol::Bool
    effectivity_steps::Int

    # Dynamics parameters
    μ::Float64
    DU::Float64
    TU::Float64
    g::Float64

    # Spacecraft parameters
    tMax::Float64
    Isp::Float64
    m0::Float64

    # Integrator settings
    maxRevs::Float64
    step::Float64
    t0::Float64
end

# Arguments in standard SI units with angles in degrees
function qLawParams(oei, oet; # [km, n.d., deg, deg, deg, deg]
    oeSwitch        = [true, true, false, false, false], 
    oeTols          = [10.0, 0.01, 0.01, 0.01, 0.01], # [km, n.d., deg, deg, deg]
    oeW             = [1.0, 1.0, 1.0, 1.0, 1.0],
    Wp              = 1.0,
    rpmin           = 6500.0,
    k               = 1000.0,
    m_petro         = 3.0,
    n_petro         = 4.0,
    r_petro         = 2.0,
    b_petro         = 0.01,
    effectivity     = false,
    abs_effectivity = false,
    rel_effectivity = false,
    η_rel_tol       = 0.1,
    η_abs_tol       = 0.1,
    effectivity_steps = 60,
    μ               = 3.9860047e5,
    DU              = 1.0,
    TU              = 1.0,
    g               = 9.80665,
    tMax            = 1.0,
    Isp             = 2000.0,
    m0              = 2000.0,
    maxRevs         = 200.0,
    step            = 5.0,
    t0              = 0.0)

    # Define constants
    deg2rad     = pi / 180.0

    # Check args
    if length(oei) != 6
        throw(ArgumentError("Vector of initial orbital elements must be of length 6."))
    end
    if length(oet) != 5
        throw(ArgumentError("Vector of target orbital elements must be of length 5."))
    end
    if length(oeSwitch) != 5
        throw(ArgumentError("Vector of target orbital element switches must be of length 5."))
    end
    if length(oeTols) != 5
        throw(ArgumentError("Vector of orbital element tolderances must be of length 5."))
    end
    if length(oeW) != 5
        throw(ArgumentError("Vector of orbital element weights must be of length 5."))
    end

    # Scale to specified units
    oei[1]      /= DU
    oei[3:6]    .*= deg2rad
    oet[1]      /= DU 
    oet[3:5]    .*= deg2rad
    oeTols[1]   /= DU 
    oeTols[3:5] .*= deg2rad
    rpmin       /= DU
    μ           *= (TU*TU / (DU*DU*DU))
    g           *= (TU*TU / (1000.0 * DU)) 
    tMax        *= (TU*TU / (1000.0 * DU))
    Isp         /= TU
    step        *= deg2rad
    t0          /= TU

    # Construct object
    qLawParams(oei,oet,oets,oetols,oeW, Wp, rpmin,k,
        m_petro,n_petro,r_petro,b_petro,effectivity,
        abs_effectivity,rel_effectivity,η_rel_tol,η_abs_tol,
        effectivity_steps,μ,DU,TU,g,tMax,Isp,m0,maxRevs,
        step, t0)
end