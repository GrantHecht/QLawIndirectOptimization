# Define some qlaw parameters
mutable struct qLawParams
    oet::Vector{Float64}
    oeW::Vector{Float64}
    Wp::Float64
    rpmin::Float64
    k::Float64
    μ::Float64
    tMax::Float64
    ηr::Float64
    steps::Int
end

function qLawParams(oet,oeW,μ,tMax)
    Wp      = 0.0
    k       = 0.0
    rpmin   = 0.0
    ηr      = 0.435
    steps   = 60
    qLawParams(oet,oeW,Wp,rpmin,k,μ,tMax,ηr,steps)
end

function qLaw(u,ps::qLawParams)
    # Compute the optimal thrust direction angles
    as = qLawThrustAngles(u[1],u[2],u[3],u[5],u[4],u[6],
            ps.oet[1],ps.oet[2],ps.oet[3],ps.oet[4],ps.oet[5],
            ps.μ,ps.tMax,ps.Wp,ps.oeW[1],ps.oeW[2],
            ps.oeW[3],ps.oeW[5],ps.oeW[4],ps.rpmin,ps.k)

    # Compute max and min Qn
    θs = range(0.0, 2*pi, length = ps.steps)
    Qnmin = Inf
    Qnmax = -Inf
    @inbounds for i in eachindex(θs)
        dQ = Qn(u[1],u[2],u[3],u[5],u[4],θs[i],
                ps.oet[1],ps.oet[2],ps.oet[3],ps.oet[5],ps.oet[4],
                ps.μ,ps.tMax,as[1],as[2],ps.Wp,ps.oeW[1],ps.oeW[2],ps.oeW[3],
                ps.oeW[5],ps.oeW[4],ps.rpmin,ps.k)
        if dQ < Qnmin; Qnmin = dQ; end
        if dQ > Qnmax; Qnmax = dQ; end
    end
    dQ = Qn(u[1],u[2],u[3],u[5],u[4],u[6],
            ps.oet[1],ps.oet[2],ps.oet[3],ps.oet[5],ps.oet[4],
            ps.μ,ps.tMax,as[1],as[2],ps.Wp,ps.oeW[1],ps.oeW[2],ps.oeW[3],
            ps.oeW[5],ps.oeW[4],ps.rpmin,ps.k)


    # Compute effectivity
    ηr  = (dQ - Qnmax) / (Qnmin - Qnmax) 
    
    # Compute coasting flag
    coast = ηr < ps.ηr

    return (as[1],as[2],coast)
end