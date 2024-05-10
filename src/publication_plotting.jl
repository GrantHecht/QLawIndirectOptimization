# Global variables to simplicy plotting for publication
const global pt_per_inch = 72
const global rad_to_deg = 180.0 / pi

# Publication figure generation
function generate_publication_figures(
    dir_path, cache::QLawTransferCache, ps::qLawParams;
)
    generate_thrust_angle_broken_axis_figure(dir_path, cache, ps)
    generate_keplerian_element_figure(dir_path, cache, ps)
    return nothing
end

function generate_thrust_angle_broken_axis_figure(dir_path, cache::QLawTransferCache, ps::qLawParams)
    CM.with_theme(CM.theme_latexfonts()) do
        # Get min and max sun angles
        not_nan = @. !isnan(cache.sun_angles)
        ymin = max(minimum(cache.sun_angles[not_nan]) - 5.0, 0.0)
        ymax = min(maximum(cache.sun_angles[not_nan]) + 5.0, 180.0)

        fig = CM.Figure(; 
            size = (3*pt_per_inch, 3*pt_per_inch), 
            fontsize = 8, 
            figure_padding = (2,10,4,4),
        )

        lims = CM.Observable(((ymin, 52.5), (127.5, ymax)))
        g = fig[1, 1] = CM.GridLayout()

        # Define base axis for y-axis label 
        CM.Axis(
            fig[1,1]; ylabel = L"$\psi$, deg",
            ytickcolor = :white, yticklabelcolor = :white, ygridvisible = false,
            xtickcolor = :white, xticklabelcolor = :white, xgridvisible = false,
            leftspinecolor = :white, rightspinecolor = :white,
            bottomspinecolor = :white, topspinecolor   = :white,
        ) 

        # Axes used for data
        ax_top = CM.Axis(fig[1, 1][1, 1])
        ax_bottom = CM.Axis(fig[1, 1][2, 1]; xlabel = "time, days")

        CM.on(lims) do (bottom, top)
            CM.ylims!(ax_bottom, bottom)
            CM.ylims!(ax_top, top)
            CM.rowsize!(g, 1, CM.Auto(top[2] - top[1]))
            CM.rowsize!(g, 2, CM.Auto(bottom[2] - bottom[1]))
        end

        CM.hidexdecorations!(ax_top, grid = false)
        ax_top.bottomspinevisible = false
        ax_bottom.topspinevisible = false

        CM.linkxaxes!(ax_top, ax_bottom)
        CM.rowgap!(g, 5)

        angle = pi/8
        linelength = 10

        segments = CM.lift(
                CM.@lift($(ax_top.yaxis.attributes.endpoints)[1]),
                CM.@lift($(ax_bottom.yaxis.attributes.endpoints)[2]),
                CM.@lift($(ax_top.elements[:yoppositeline][1])[1]),
                CM.@lift($(ax_bottom.elements[:yoppositeline][1])[2]),
            ) do p1, p2, p3, p4
            ps = CM.Point2f[p1, p2, p3, p4]
            
            map(ps) do p
                a = p + CM.Point2f(cos(angle), sin(angle)) * 0.5 * linelength
                b = p - CM.Point2f(cos(angle), sin(angle)) * 0.5 * linelength
                (a, b)
            end
        end

        # segments = CM.lift(
        #         CM.@lift($(ax_top.yaxis.attributes.endpoints)[1]),
        #         CM.@lift($(ax_bottom.yaxis.attributes.endpoints)[2]),
        #     ) do p1, p2
        #     ps = CM.Point2f[p1, p2]
            
        #     map(ps) do p
        #         a = p + CM.Point2f(cos(angle), sin(angle)) * 0.5 * linelength
        #         b = p - CM.Point2f(cos(angle), sin(angle)) * 0.5 * linelength
        #         (a, b)
        #     end
        # end

        CM.linesegments!(fig.scene, segments; color = :black, linewidth = 0.75)

        # Set ticks
        ax_bottom.yticks = 0:10:50
        ax_top.yticks = 130:10:180

        # Strip data to reduce points plotted
        times_stripped = cache.times[1:10:end]
        angles_stripped = cache.sun_angles[1:10:end]

        below_60 = @. angles_stripped < 60.0
        above_120 = @. angles_stripped >= 120.0

        CM.scatter!(
            ax_bottom, times_stripped[below_60], angles_stripped[below_60];
            markersize = 1,
        )
        CM.scatter!(
            ax_top, times_stripped[above_120], angles_stripped[above_120];
            markersize = 1,
        )
        CM.notify(lims)

        # Save figure
        CM.save(joinpath(dir_path, "sun_angles.png"), fig)
    end
    return nothing
end

function generate_thrust_angle_figure(dir_path, cache::QLawTransferCache, ps::qLawParams)
    CM.with_theme(CM.theme_latexfonts()) do
        # Instantiate figure
        fig = CM.Figure(; size = (4*pt_per_inch, 3*pt_per_inch), fontsize = 8)

        # Define main axis
        ax_main = CM.Axis(
            fig[1,1]; 
            xlabel = "time, days", 
            ylabel = L"$\psi$, deg",
        )

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
    # Get distance unit
    DU = ps.meePs.LU

    # Convert cartesian states to keplerian elements
    kep_states = zeros(length(cache.states), 7)
    for i in eachindex(cache.states)
        kep = AstroUtils.convertState(cache.states[i], AstroUtils.Cartesian, AstroUtils.Keplerian, ps.μ)[1]
        kep_states[i, 1:6] .= kep
        kep_states[i, 7] = cache.states[i][7]

        # Convert to degrees
        kep_states[i, 1] *= DU
        kep_states[i, 3:6] .*= rad_to_deg
    end

    # Define verticle size of each exis
    ax_vert_size = 1.1

    # Determine which states we need to plot
    p_flags = ps.oeW .> 0.0
    n_plots = sum(p_flags) + 1

    CM.with_theme(CM.theme_latexfonts()) do
        # Instantiate figure
        fig = CM.Figure(; size = (4*pt_per_inch, n_plots*ax_vert_size*pt_per_inch), fontsize = 8, figure_padding = 1)

        ax_kwargs = (;)

        ax_count = 0
        axs = []
        if p_flags[1]
            ax_count += 1
            # Semi-major axis
            ax_sma = CM.Axis(fig[ax_count,1]; ylabel = L"$a$, km", ax_kwargs...)
            CM.hlines!(ax_sma, [DU*ps.oet[1]]; color = :black, linewidth = 0.5)
            CM.lines!(ax_sma, cache.times, kep_states[:, 1]; linewidth = 0.5)
            push!(axs, ax_sma)
        end

        if p_flags[2]
            ax_count += 1
            # eccentricity axis
            ax_e = CM.Axis(fig[ax_count,1]; ylabel = L"$e$", ax_kwargs...)
            CM.hlines!(ax_e, [ps.oet[2]]; color = :black, linewidth = 0.5)
            CM.lines!(ax_e, cache.times, kep_states[:, 2]; linewidth = 0.5)
            push!(axs, ax_e)
        end

        if p_flags[3]
            ax_count += 1
            # inclination axis
            ax_inc = CM.Axis(fig[ax_count,1]; ylabel = L"$i$, deg", ax_kwargs...)
            CM.hlines!(ax_inc, [rad_to_deg*ps.oet[3]]; color = :black, linewidth = 0.5)
            CM.lines!(ax_inc, cache.times, kep_states[:, 3]; linewidth = 0.5)
            push!(axs, ax_inc)
        end

        if p_flags[4]
            ax_count += 1
            # RAAN
            ax_raan = CM.Axis(fig[ax_count,1]; ylabel = L"$\Omega$, deg", ax_kwargs...)
            CM.hlines!(ax_raan, [rad_to_deg*ps.oet[4]]; color = :black, linewidth = 0.5)
            CM.lines!(ax_raan, cache.times, kep_states[:, 4]; linewidth = 0.5)
            push!(axs, ax_raan)
        end

        if p_flags[5]
            ax_count += 1
            # arg. of pariapsis axis
            ax_aop = CM.Axis(fig[ax_count,1]; ylabel = L"$\omega$, deg", ax_kwargs...)
            CM.hlines!(ax_aop, [rad_to_deg*ps.oet[5]]; color = :black, linewidth = 0.5)
            CM.lines!(ax_aop, cache.times, kep_states[:, 5]; linewidth = 0.5)
            push!(axs, ax_aop)
        end

        # mass axis
        ax_count += 1
        ax_m = CM.Axis(fig[ax_count,1]; xlabel = L"time, days$$", ylabel = L"$m$, kg", ax_kwargs...)
        CM.lines!(ax_m, cache.times, kep_states[:, 7]; linewidth = 0.5)
        push!(axs, ax_m)

        # Link axes
        CM.linkxaxes!(axs...)

        # Save figure
        CM.save(joinpath(dir_path, "keplerian_elements.png"), fig)
    end
    return nothing
end