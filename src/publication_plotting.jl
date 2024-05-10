# Global variables to simplicy plotting for publication
const global pt_per_inch = 72
const global rad_to_deg = 180.0 / pi

# Publication figure generation
function generate_publication_figures(
    dir_path, cache::QLawTransferCache, ps::qLawParams;
)
    generate_thrust_angle_figure(dir_path, cache, ps)
    generate_keplerian_element_figure(dir_path, cache, ps)
    return nothing
end

function generate_thrust_angle_figure(dir_path, cache::QLawTransferCache, ps::qLawParams)
    if !all(isnan.(cache.sun_angles))
        # Instantiate figure
        fig = CM.Figure(; size = (4*pt_per_inch, 3*pt_per_inch), fontsize = 10)
        ax = CM.Axis(fig[1,1]; xlabel = L"time, days$$", ylabel = L"$\psi$, deg")

        # Strip data to reduce points plotted
        times_stripped = cache.times[1:10:end]
        angles_stripped = cache.sun_angles[1:10:end]

        # Plot sun-angles
        CM.scatter!(
            ax, times_stripped, angles_stripped;
            markersize = 2,
        )

        # Save figure
        CM.save(joinpath(dir_path, "sun_angles.png"), fig)
    end
    return nothing
end

function generate_keplerian_element_figure(dir_path, cache::QLawTransferCache, ps::qLawParams)
    # Convert cartesian states to keplerian elements
    kep_states = zeros(length(cache.states), 7)
    for i in eachindex(cache.states)
        kep = AstroUtils.convertState(cache.states[i], AstroUtils.Cartesian, AstroUtils.Keplerian, ps.μ)[1]
        kep_states[i, 1:6] .= kep
        kep_states[i, 7] = cache.states[i][7]

        # Convert to degrees
        kep_states[i, 3:6] .*= rad_to_deg
    end

    # Define verticle size of each exis
    ax_vert_size = 1.25

    # Determine which states we need to plot
    p_flags = ps.oeW .> 0.0
    n_plots = sum(p_flags) + 1

    # Instantiate figure
    fig = CM.Figure(; size = (4*pt_per_inch, n_plots*ax_vert_size*pt_per_inch), fontsize = 10)

    # Axis kwargs
    ax_kwargs = (; yminorticksvisible = true, yminorgridvisible = true)

    ax_count = 0
    if p_flags[1]
        ax_count += 1
        # Semi-major axis
        ax_sma = CM.Axis(fig[ax_count,1]; yminorticks = CM.IntervalsBetween(2), xlabel = L"time, days$$", ylabel = L"$a$, km", ax_kwargs...)
        CM.lines!(ax_sma, cache.times, kep_states[:, 1])
        CM.hlines!(ax_sma, [ps.oet[1]])
    end

    if p_flags[2]
        ax_count += 1
        # eccentricity axis
        ax_e = CM.Axis(fig[ax_count,1]; yminorticks = CM.IntervalsBetween(2), xlabel = L"time, days$$", ylabel = L"$e$", ax_kwargs...)
        CM.lines!(ax_e, cache.times, kep_states[:, 2])
        CM.hlines!(ax_e, [ps.oet[2]])
    end

    if p_flags[3]
        ax_count += 1
        # inclination axis
        ax_inc = CM.Axis(fig[ax_count,1]; yminorticks = CM.IntervalsBetween(2), xlabel = L"time, days$$", ylabel = L"$i$, deg", ax_kwargs...)
        CM.lines!(ax_inc, cache.times, kep_states[:, 3])
        CM.hlines!(ax_inc, [rad_to_deg*ps.oet[3]])
    end

    if p_flags[4]
        ax_count += 1
        # RAAN
        ax_inc = CM.Axis(fig[ax_count,1]; yminorticks = CM.IntervalsBetween(2), xlabel = L"time, days$$", ylabel = L"$\Omega$, deg", ax_kwargs...)
        CM.lines!(ax_inc, cache.times, kep_states[:, 4])
        CM.hlines!(ax_inc, [rad_to_deg*ps.oet[4]])
    end

    if p_flags[5]
        ax_count += 1
        # arg. of pariapsis axis
        ax_inc = CM.Axis(fig[ax_count,1]; yminorticks = CM.IntervalsBetween(2), xlabel = L"time, days$$", ylabel = L"$\omega$, deg", ax_kwargs...)
        CM.lines!(ax_inc, cache.times, kep_states[:, 5])
        CM.hlines!(ax_inc, [rad_to_deg*ps.oet[5]])
    end

    # mass axis
    ax_count += 1
    ax_m = CM.Axis(fig[ax_count,1]; yminorticks = CM.IntervalsBetween(2), xlabel = L"time, days$$", ylabel = L"$m$, kg", ax_kwargs...)
    CM.lines!(ax_m, cache.times, kep_states[:, 7])

    # Save figure
    CM.save(joinpath(dir_path, "keplerian_elements.png"), fig)
    return nothing
end