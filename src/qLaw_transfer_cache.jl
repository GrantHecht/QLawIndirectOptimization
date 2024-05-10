
struct QLawTransferCache
    # Fields for full transfer trajectory
    times::Vector{Float64}
    states::Vector{SVector{7,Float64}}
    angles::Vector{SVector{2,Float64}}
    thrust_vals::Vector{Float64}
    sun_angles::Vector{Float64}
    coast_positive_dQs::Vector{Bool}
    coast_effectives::Vector{Bool}
    coast_eclipse::Vector{Bool}
    thrid_body_perturbations::Vector{SVector{3,Float64}}
    J2_perturbations::Vector{SVector{3,Float64}}

    # Fields for state at steps
    time_at_steps::Vector{Float64}
    state_at_steps::Vector{SVector{7,Float64}}

    # Fields for qLaw parameters
    weights::Vector{Float64}
    effectivity_tol::Float64

    function QLawTransferCache(weights::Vector{Float64}, η::Float64)
        return new(
            Vector{Float64}(undef, 0),
            Vector{SVector{6,Float64}}(undef, 0),
            Vector{SVector{2,Float64}}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Bool}(undef, 0),
            Vector{Bool}(undef, 0),
            Vector{Bool}(undef, 0),
            Vector{SVector{3,Float64}}(undef, 0),
            Vector{SVector{3,Float64}}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{SVector{7,Float64}}(undef, 0),
            weights, η,
        )
    end
end

function dump_to_mat(cache::QLawTransferCache, filename::String)
    # Construct matricies from Vector{SVector{...}}(s)
    states_mat = zeros(length(cache.states), 7)
    angles_mat = zeros(length(cache.angles), 2)
    states_at_steps_mat = zeros(length(cache.state_at_steps), 7)
    for i in eachindex(cache.states)
        states_mat[i,:] .= cache.states[i]
        angles_mat[i,:] .= cache.angles[i]
    end
    tbody_mat  = zeros(length(cache.thrid_body_perturbations), 3)
    for i in eachindex(cache.thrid_body_perturbations)
        tbody_mat[i,:] .= cache.thrid_body_perturbations[i]
    end
    J2_mat     = zeros(length(cache.J2_perturbations), 3)
    for i in eachindex(cache.J2_perturbations)
        J2_mat[i,:] .= cache.J2_perturbations[i]
    end
    for i in eachindex(cache.state_at_steps)
        states_at_steps_mat[i,:] .= cache.state_at_steps[i]
    end

    # Create dictionary
    dict = Dict(
        "times" => cache.times,
        "states" => states_mat,
        "angles" => angles_mat,
        "thrust_vals" => cache.thrust_vals,
        "sun_angles" => cache.sun_angles,
        "coast_positive_dQs" => cache.coast_positive_dQs,
        "coast_effectives" => cache.coast_effectives,
        "coast_eclipse" => cache.coast_eclipse,
        "thrid_body_perturbations" => tbody_mat,
        "J2_perturbations" => J2_mat,
        "time_at_steps" => cache.time_at_steps,
        "state_at_steps" => states_at_steps_mat,
        "weights" => cache.weights,
        "effectivity_tol" => cache.effectivity_tol,
    )

    # Write to file
    MAT.matwrite(filename, dict)

    return nothing
end

function push_update!(cache::QLawTransferCache, de_sol, α, β, thrust, eclipse_coast, effectivity_coast, positive_dQ_coast, ps::qLawParams)
    # Check if cache is empty
    if length(cache.times) == 0
        state_mee_s = de_sol.u[1]
        time        = state_mee_s[6]
        state_mee   = SA[state_mee_s[1], state_mee_s[2], state_mee_s[3], state_mee_s[4], state_mee_s[5], de_sol.t[1]]
        state_cart  = AstroUtils.convertState(state_mee, AstroUtils.MEE, AstroUtils.Cartesian, ps.μ)
        state       = SA[state_cart[1], state_cart[2], state_cart[3], state_cart[4], state_cart[5], state_cart[6], state_mee_s[7]]
        push!(cache.time_at_steps,  time)
        push!(cache.state_at_steps, state) 
    end

    # Update cache from differential equation solution
    Ls = range(de_sol.t[1], de_sol.t[end]; length = ps.savedStatesAtSteps)
    for L in Ls
        # Get state
        state_mee_s = de_sol(L)
        time        = state_mee_s[6]
        state_mee   = SA[state_mee_s[1], state_mee_s[2], state_mee_s[3], state_mee_s[4], state_mee_s[5], L]
        state_cart  = AstroUtils.convertState(state_mee, AstroUtils.MEE, AstroUtils.Cartesian, ps.μ)
        state       = SA[state_cart[1], state_cart[2], state_cart[3], state_cart[4], state_cart[5], state_cart[6], state_mee_s[7]]

        # Push trivial info
        push!(cache.times,  time)
        push!(cache.states, state)
        push!(cache.angles, SA[α, β])
        push!(cache.thrust_vals, thrust)
        push!(cache.coast_eclipse,      eclipse_coast)
        push!(cache.coast_effectives,   effectivity_coast)
        push!(cache.coast_positive_dQs, positive_dQ_coast)

        # Push perturbations
        if ps.savePerturbations
            if ps.meePs.thirdBodyPerterbations
                tb_pert = get_third_body_perturbation(state_cart, time, ps)
                push!(cache.thrid_body_perturbations, tb_pert)
            end
            if ps.meePs.onlyJ2
                j2_pert = get_j2_perturbation(state_mee, ps)
                push!(cache.J2_perturbations, j2_pert)
            end
        end

        # Compute sun-angle
        if ps.meePs.thirdBodyEphemerides !== nothing && 10 in ps.meePs.thirdBodyEphemerides.targIDs
            # Position vectors
            rs  = AstroUtils.getPosition(ps.meePs.thirdBodyEphemerides, 10, ps.meePs.initEpoch + time*ps.meePs.TU)
            rss = SA[
                rs[1]/ps.meePs.LU - state[1], 
                rs[2]/ps.meePs.LU - state[2], 
                rs[3]/ps.meePs.LU - state[3],
            ]

            # Rotation matrix
            RlI = get_lvlh_to_inertial_mee(state_mee, ps)
            
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


function plot_transfer(
    cache::QLawTransferCache, ps::qLawParams;
    axes = SA[1,2],
    show_coast = true,
    linewidth = 1.0,
    linewidth_terminal = 1.0,
)
    if show_coast
        return plot_transfer_with_coasts(cache, ps, axes; linewidth = linewidth, linewidth_terminal = linewidth_terminal)
    else
        return plot_transfer_without_coasts(cache, ps, axes)
    end
end

function plot_transfer(
    file::String, cache::QLawTransferCache, ps::qLawParams;
    axes = SA[1,2],
    show_coast = true,
    linewidth = 1.0,
    linewidth_terminal = 1.0,
)
    fig = plot_transfer(
        cache, ps; 
        axes = axes, 
        show_coast = show_coast, 
        linewidth = linewidth,
        linewidth_terminal = linewidth_terminal,
    )
    CM.save(file, fig)
    return nothing
end


function plot_transfer_with_coasts(
    cache::QLawTransferCache, ps::qLawParams, axes;
    linewidth = 1.0, 
    linewidth_terminal = 1.0,
)
    # Grab length unit
    DU = ps.meePs.LU

    # Grab the thrust and coast arcs
    a1_ta, a2_ta, a1_ef, a2_ef, a1_pq, a2_pq, a1_ec, a2_ec = get_thrust_and_coast_arcs(cache, ps, axes)

    # Get initial and final orbit
    x0 = cache.states[1]
    xf = cache.states[end]
    a1_ik, a2_ik = get_keplerian_orbit(x0, ps, axes)
    a1_fk, a2_fk = get_keplerian_orbit(xf, ps, axes)

    # Plot
    fig = CM.Figure(); ax = CM.Axis(fig[1,1]; aspect=CM.DataAspect())
    CM.lines!(ax, a1_ta, a2_ta, label="Thrust Arcs", color=:red, linewidth=linewidth)
    CM.lines!(ax, a1_ec, a2_ec, label="Ecc. Coast Arcs", color=:green, linewidth=linewidth)
    CM.lines!(ax, a1_ef, a2_ef, label="Eff. Coast Arcs", color=:blue, linewidth=linewidth)
    CM.lines!(ax, a1_pq, a2_pq, label="dQ > 0 Coast Arcs", color=:deepskyblue, linewidth=linewidth)
    CM.lines!(ax, a1_ik, a2_ik, label="Initial Orbit", color=:black, linewidth=linewidth_terminal)
    CM.lines!(ax, a1_fk, a2_fk, label="Final Orbit", color=:black, linewidth=linewidth_terminal)
    return fig
end

function plot_transfer_without_coasts(
    cache::QLawTransferCache, ps::qLawParams, axes
)

end

function get_keplerian_orbit(cart, ps::qLawParams, axes; n::Int = 1000)
    # Get parameters
    μ = ps.μ
    DU = ps.meePs.LU

    # Handle states
    r0 = SA[cart[1], cart[2], cart[3]]
    v0 = SA[cart[4], cart[5], cart[6]]

    # Compute orital period
    ξ = 0.5*norm(v0)^2 - μ/norm(r0)
    a = -μ / (2.0*ξ)
    P = 2.0*π*sqrt(a^3 / μ)

    # Compute orbit
    Δts = range(0.0, P; length = n)
    first_axis_pos = Vector{Float32}(undef, n)
    second_axis_pos = Vector{Float32}(undef, n)
    @inbounds for i in eachindex(Δts)
        Δt = Δts[i]
        r, v, flag = kepler(r0, v0, Δt, μ)
        first_axis_pos[i] = DU*r[axes[1]]
        second_axis_pos[i] = DU*r[axes[2]]
    end
    return first_axis_pos, second_axis_pos
end

function get_thrust_and_coast_arcs(
    cache::QLawTransferCache, ps::qLawParams, axes
)
    # Get distance unit
    DU = ps.meePs.LU

    n = length(cache.times)
    first_axis_thrust_arcs  = Vector{Float32}(undef, n)
    second_axis_thrust_arcs  = Vector{Float32}(undef, n)
    first_axis_eff_coast_arcs = Vector{Float32}(undef, n)
    second_axis_eff_coast_arcs = Vector{Float32}(undef, n)
    first_axis_posqd_coast_arcs = Vector{Float32}(undef, n)
    second_axis_posqd_coast_arcs = Vector{Float32}(undef, n)
    first_axis_eclipse_coast_arcs = Vector{Float32}(undef, n)
    second_axis_eclipse_coast_arcs = Vector{Float32}(undef, n)
    @inbounds for i in eachindex(cache.times)
        coast = cache.coast_eclipse[i] || cache.coast_effectives[i] || cache.coast_positive_dQs[i]
        if !coast
            first_axis_thrust_arcs[i] = DU*cache.states[i][axes[1]]
            second_axis_thrust_arcs[i] = DU*cache.states[i][axes[2]]
            first_axis_eff_coast_arcs[i] = NaN
            second_axis_eff_coast_arcs[i] = NaN
            first_axis_posqd_coast_arcs[i] = NaN
            second_axis_posqd_coast_arcs[i] = NaN
            first_axis_eclipse_coast_arcs[i] = NaN
            second_axis_eclipse_coast_arcs[i] = NaN
        else
            if cache.coast_eclipse[i]
                first_axis_thrust_arcs[i] = NaN
                second_axis_thrust_arcs[i] = NaN
                first_axis_eff_coast_arcs[i] = NaN
                second_axis_eff_coast_arcs[i] = NaN
                first_axis_posqd_coast_arcs[i] = NaN
                second_axis_posqd_coast_arcs[i] = NaN
                first_axis_eclipse_coast_arcs[i] = DU*cache.states[i][axes[1]]
                second_axis_eclipse_coast_arcs[i] = DU*cache.states[i][axes[2]]
            elseif cache.coast_effectives[i]
                first_axis_thrust_arcs[i] = NaN
                second_axis_thrust_arcs[i] = NaN
                first_axis_eff_coast_arcs[i] = DU*cache.states[i][axes[1]]
                second_axis_eff_coast_arcs[i] = DU*cache.states[i][axes[2]]
                first_axis_posqd_coast_arcs[i] = NaN
                second_axis_posqd_coast_arcs[i] = NaN
                first_axis_eclipse_coast_arcs[i] = NaN
                second_axis_eclipse_coast_arcs[i] = NaN
            else
                first_axis_thrust_arcs[i] = NaN
                second_axis_thrust_arcs[i] = NaN
                first_axis_eff_coast_arcs[i] = NaN
                second_axis_eff_coast_arcs[i] = NaN
                first_axis_posqd_coast_arcs[i] = DU*cache.states[i][axes[1]]
                second_axis_posqd_coast_arcs[i] = DU*cache.states[i][axes[2]]
                first_axis_eclipse_coast_arcs[i] = NaN
                second_axis_eclipse_coast_arcs[i] = NaN
            end
        end
    end
    return (
        first_axis_thrust_arcs,
        second_axis_thrust_arcs,
        first_axis_eff_coast_arcs,
        second_axis_eff_coast_arcs,
        first_axis_posqd_coast_arcs,
        second_axis_posqd_coast_arcs,
        first_axis_eclipse_coast_arcs,
        second_axis_eclipse_coast_arcs,
    )
end



