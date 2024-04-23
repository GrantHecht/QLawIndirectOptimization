
struct QLawTransferCache
    # Fields for full transfer trajectory
    times::Vector{Float64}
    states::Vector{SVector{7,Float64}}
    angles::Vector{SVector{2,Float64}}
    thrust_vals::Vector{Float64}
    sun_angles::Vector{Float64}
    coast_positive_dQs::Vector{Bool}
    coast_effectives::Vector{Bool}

    # Fields for state at steps
    time_at_steps::Vector{Float64}
    state_at_steps::Vector{SVector{7,Float64}}

    function QLawTransferCache()
        return new(
            Vector{Float64}(undef, 0),
            Vector{SVector{6,Float64}}(undef, 0),
            Vector{SVector{2,Float64}}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Bool}(undef, 0),
            Vector{Bool}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{SVector{7,Float64}}(undef, 0),
        )
    end
end

function push_update!(cache::QLawTransferCache, de_sol, α, β, thrust, effectivity_coast, positive_dQ_coast, ps::qLawParams)
    # Check if cache is empty
    if length(cache.times) == 0
        state_sund = de_sol.u[1]
        time  = state_sund[6]
        state = SA[state_sund[1], state_sund[2], state_sund[3], state_sund[4], state_sund[5], de_sol.t[1], state_sund[7]]
        push!(cache.time_at_steps,  time)
        push!(cache.state_at_steps, state) 
    end

    # Update cache from differential equation solution
    for i in eachindex(de_sol)
        # Get state
        state_sund = de_sol.u[i]
        time  = state_sund[6]
        state = SA[state_sund[1], state_sund[2], state_sund[3], state_sund[4], state_sund[5], de_sol.t[i], state_sund[7]]

        # Push trivial info
        push!(cache.times,  time)
        push!(cache.states, state)
        push!(cache.angles, SA[α, β])
        push!(cache.thrust_vals, thrust)
        push!(cache.coast_effectives,   effectivity_coast)
        push!(cache.coast_positive_dQs, positive_dQ_coast)

        # Compute sun-angle
        if ps.meePs.thirdBodyEphemerides !== nothing && 10 in ps.meePs.thirdBodyEphemerides.targIDs
            # Position vectors
            rs  = AstroUtils.getPosition(ps.meePs.thirdBodyEphemerides, 10, ps.meePs.initEpoch + time*ps.meePs.TU)
            rsc = AstroUtils.convertState(state[SA[1,2,3,4,5,6]], AstroUtils.MEE, AstroUtils.Cartesian, ps.μ)
            rss = SA[rs[1]/ps.meePs.LU - rsc[1], rs[2]/ps.meePs.LU - rsc[2], rs[3]/ps.meePs.LU - rsc[3]]

            # Rotation matrix
            RlI = get_lvlh_to_inertial(state[SA[1,2,3,4,5,6]], ps)
            
            # Control vector
            ut  = RlI*SA[cos(β)*sin(α), cos(β)*cos(α), sin(β)]

            # Push sun angle
            push!(cache.sun_angles, acosd(dot(ut, rss) / (norm(ut)*norm(rss))))
        else
            push!(cache.sun_angles, NaN)
        end
    end

    # Update states at steps
    push!(cache.time_at_steps,  de_sol.t[end])
    push!(cache.state_at_steps, de_sol.u[end])
end

function plot_transfer(cache::QLawTransferCache, ps::qLawParams)
    fig = Figure();
    ax = Axis3(fig[1, 1], aspect = :data)
    xs = zeros(length(cache.times))
    ys = zeros(length(cache.times))
    zs = zeros(length(cache.times))
    for i in eachindex(cache.times)
        cart = AstroUtils.convertState(cache.states[i], AstroUtils.MEE, AstroUtils.Cartesian, ps.μ) 
        xs[i] = cart[1]
        ys[i] = cart[2]
        zs[i] = cart[3]
    end
    lines!(ax, xs, ys, zs)
    return fig
end