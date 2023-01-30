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
    α::Float64
    β::Float64
    coasting::Bool
end

function qLawParams(oet,oeW,Wp,rpmin,k,μ,tMax,ηr,steps)
    return qLawParams(oet,oeW,Wp,rpmin,k,μ,tMax,ηr,steps,0.0,0.0,false)
end