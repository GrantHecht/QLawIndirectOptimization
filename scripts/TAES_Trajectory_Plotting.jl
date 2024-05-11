
using QLawIndirectOptimization
using JLD2, FileIO
using CairoMakie
using StaticArrays

const global QP = QLawIndirectOptimization

function generate_gto_to_geo_figures(
    dir_path, cache::QP.QLawTransferCache, ps::qLawParams;
)
    # Get points per inch 
    ppi = QP.pt_per_inch

    # Define the axes we're plotting
    axes = (SA[1,2], SA[2,3])

    # Grab distance unit from parameters
    DU = ps.meePs.LU

    with_theme(theme_latexfonts()) do
        # Grab the thrust and coast arcs for axes
        a1_ta1, a2_ta1, a1_ef1, a2_ef1, a1_pq1, a2_pq1, a1_ec1, a2_ec1 = 
            QP.get_thrust_and_coast_arcs(cache, ps, axes[1])
        # Grab the thrust and coast arcs
        a1_ta2, a2_ta2, a1_ef2, a2_ef2, a1_pq2, a2_pq2, a1_ec2, a2_ec2 = 
            QP.get_thrust_and_coast_arcs(cache, ps, axes[2])

        # Get spans for axis 1
        not_nan = @. !isnan(a1_ta1)
        xmin1   = minimum(a1_ta1[not_nan])
        xmax1   = maximum(a1_ta1[not_nan])
        xspan1  = xmax1 - xmin1
        not_nan = @. !isnan(a2_ta1)
        ymin1   = minimum(a2_ta1[not_nan])
        ymax1   = maximum(a2_ta1[not_nan])
        yspan1  = ymax1 - ymin1

        # Get spans for axis 2
        not_nan = @. !isnan(a1_ta2)
        xmin2   = minimum(a1_ta2[not_nan])
        xmax2   = maximum(a1_ta2[not_nan])
        xspan2  = xmax2 - xmin2
        not_nan = @. !isnan(a2_ta2)
        ymin2   = minimum(a2_ta2[not_nan])
        ymax2   = maximum(a2_ta2[not_nan])
        yspan2  = ymax2 - ymin2

        # Handle x limits
        xmin    = min(xmin1, xmin2) - 0.75e4
        xmax    = max(xmax1, xmax2) + 5e4

        # Create figure
        fig = Figure(; size = (4*ppi, 4*ppi), fontsize = 8, figure_padding = 1); 

        # Define axes and link
        ax1 = Axis(
            fig[1, 1], 
            aspect = DataAspect(), 
            xlabel = L"$x$, km", 
            ylabel = L"$y$, km", 
            #xticklabelrotation = pi/4,
        )  
        ax2 = Axis(
            fig[2, 1], 
            aspect = DataAspect(), 
            xlabel = L"$y$, km", 
            ylabel = L"$z$, km", 
            #xticklabelrotation = pi/4,
        )  
        linkxaxes!(ax1, ax2)

        # ==== First axis
        # Get the initial and final orbits
        x0 = cache.states[1]
        xf = cache.states[end]
        a1_ik, a2_ik = QP.get_keplerian_orbit(x0, ps, axes[1])
        a1_fk, a2_fk = QP.get_keplerian_orbit(xf, ps, axes[1])

        lines!(ax1, a1_ta1, a2_ta1, label="Thrusting", color=:red, linewidth=0.25)
        lines!(ax1, a1_ec1, a2_ec1, label="Eclipsed", color=:green, linewidth=0.25)
        if ps.ηr != 0.0
            lines!(ax1, a1_ef1, a2_ef1, label="Coasting", color=:blue, linewidth=0.25)
        end
        lines!(ax1, a1_ik, a2_ik, label="Initial Orbit", color=:black, linewidth=1.0)
        lines!(ax1, a1_fk, a2_fk, label="Target Orbit", color=:black, linewidth=1.0)

        # ==== Second axis
        # Get the initial and final orbits
        x0 = cache.states[1]
        xf = cache.states[end]
        a1_ik, a2_ik = QP.get_keplerian_orbit(x0, ps, axes[2])
        a1_fk, a2_fk = QP.get_keplerian_orbit(xf, ps, axes[2])

        lines!(ax2, a1_ta2, a2_ta2, label="Thrusting", color=:red, linewidth=0.25)
        lines!(ax2, a1_ec2, a2_ec2, label="Eclipsed", color=:green, linewidth=0.25)
        if ps.ηr != 0.0
            lines!(ax2, a1_ef2, a2_ef2, label="Coasting", color=:blue, linewidth=0.25)
        end
        lines!(ax2, a1_ik, a2_ik, label="Initial Orbit", color=:black, linewidth=1.0)
        lines!(ax2, a1_fk, a2_fk, label="Target Orbit", color=:black, linewidth=1.0)

        xlims!(ax1, xmin, xmax)
        xlims!(ax2, xmin, xmax)
        rowsize!(fig.layout, 1, Auto(yspan1))
        rowsize!(fig.layout, 2, Auto(yspan2))

        # Define legend
        if ps.ηr == 0.0
            axislegend(
                ax1, 
                [LineElement(color=:red), LineElement(color=:green)], 
                ["Thrusting", "Eclipsed"];
                position = :rt,
                height = 0.75ppi,
                framevisible = false,
            )
        else
            axislegend(
                ax1, 
                [LineElement(color=:red), LineElement(color=:blue), LineElement(color=:green)], 
                ["Thrusting", "Coasting", "Eclipsed"], 
                position = :rt,
                height = 0.8ppi,
                framevisible = false,
            )
        end

        save(joinpath(dir_path, "transfer.pdf"), fig)
    end
    return nothing
end

function generate_molniya_like_figures(
    dir_path, cache::QP.QLawTransferCache, ps::qLawParams;
)
    # Get points per inch 
    ppi = QP.pt_per_inch

    # Define the axes we're plotting
    axes = (SA[1,3], SA[2,3])

    # Grab distance unit from parameters
    DU = ps.meePs.LU

    with_theme(theme_latexfonts()) do
        # Grab the thrust and coast arcs for axes
        a1_ta1, a2_ta1, a1_ef1, a2_ef1, a1_pq1, a2_pq1, a1_ec1, a2_ec1 = 
            QP.get_thrust_and_coast_arcs(cache, ps, axes[1])
        # Grab the thrust and coast arcs
        a1_ta2, a2_ta2, a1_ef2, a2_ef2, a1_pq2, a2_pq2, a1_ec2, a2_ec2 = 
            QP.get_thrust_and_coast_arcs(cache, ps, axes[2])

        # Get spans for axis 1
        not_nan = @. !isnan(a1_ta1)
        xmin1   = minimum(a1_ta1[not_nan])
        xmax1   = maximum(a1_ta1[not_nan])
        xspan1  = xmax1 - xmin1
        not_nan = @. !isnan(a2_ta1)
        ymin1   = minimum(a2_ta1[not_nan])
        ymax1   = maximum(a2_ta1[not_nan])
        yspan1  = ymax1 - ymin1

        # Get spans for axis 2
        not_nan = @. !isnan(a1_ta2)
        xmin2   = minimum(a1_ta2[not_nan])
        xmax2   = maximum(a1_ta2[not_nan])
        xspan2  = xmax2 - xmin2
        not_nan = @. !isnan(a2_ta2)
        ymin2   = minimum(a2_ta2[not_nan])
        ymax2   = maximum(a2_ta2[not_nan])
        yspan2  = ymax2 - ymin2

        # Handle x limits
        xmin    = min(xmin1, xmin2) - 0.75e4
        xmax    = max(xmax1, xmax2) + 5e4

        # Handle y limits
        ymin    = min(ymin1, ymin2)
        ymax    = max(ymax2, ymax2)

        # Create figure
        fig = Figure(; size = (6*ppi, 3*ppi), fontsize = 8, figure_padding = 1); 

        # Define axes and link
        ax1 = Axis(
            fig[1, 1], 
            aspect = DataAspect(), 
            xlabel = L"$x$, km", 
            ylabel = L"$z$, km", 
        )  
        ax2 = Axis(
            fig[1, 2];
            aspect = DataAspect(), 
            xlabel = L"$y$, km", 
            ytickcolor = :white,
            yticklabelcolor = :white
        )  
        linkyaxes!(ax1, ax2)

        # ==== First axis
        # Get the initial and final orbits
        x0 = cache.states[1]
        xf = cache.states[end]
        a1_ik, a2_ik = QP.get_keplerian_orbit(x0, ps, axes[1])
        a1_fk, a2_fk = QP.get_keplerian_orbit(xf, ps, axes[1])

        lines!(ax1, a1_ta1, a2_ta1, label="Thrusting", color=:red, linewidth=0.25)
        lines!(ax1, a1_ec1, a2_ec1, label="Eclipsed", color=:green, linewidth=0.25)
        if ps.ηr != 0.0
            lines!(ax1, a1_ef1, a2_ef1, label="Coasting", color=:blue, linewidth=0.25)
        end
        lines!(ax1, a1_ik, a2_ik, label="Initial Orbit", color=:black, linewidth=1.0)
        lines!(ax1, a1_fk, a2_fk, label="Target Orbit", color=:black, linewidth=1.0)

        # ==== Second axis
        # Get the initial and final orbits
        x0 = cache.states[1]
        xf = cache.states[end]
        a1_ik, a2_ik = QP.get_keplerian_orbit(x0, ps, axes[2])
        a1_fk, a2_fk = QP.get_keplerian_orbit(xf, ps, axes[2])

        lines!(ax2, a1_ta2, a2_ta2, label="Thrusting", color=:red, linewidth=0.25)
        lines!(ax2, a1_ec2, a2_ec2, label="Eclipsed", color=:green, linewidth=0.25)
        if ps.ηr != 0.0
            lines!(ax2, a1_ef2, a2_ef2, label="Coasting", color=:blue, linewidth=0.25)
        end
        lines!(ax2, a1_ik, a2_ik, label="Initial Orbit", color=:black, linewidth=1.0)
        lines!(ax2, a1_fk, a2_fk, label="Target Orbit", color=:black, linewidth=1.0)

        ylims!(ax1, ymin, ymax)
        ylims!(ax2, ymin, ymax)
        colsize!(fig.layout, 1, Auto(xspan1))
        colsize!(fig.layout, 2, Auto(xspan2))

        # # Define legend
        # if ps.ηr == 0.0
        #     axislegend(
        #         ax1, 
        #         [LineElement(color=:red), LineElement(color=:green)], 
        #         ["Thrusting", "Eclipsed"];
        #         position = :rt,
        #         height = 0.75ppi,
        #         framevisible = false,
        #     )
        # else
        #     axislegend(
        #         ax1, 
        #         [LineElement(color=:red), LineElement(color=:blue), LineElement(color=:green)], 
        #         ["Thrusting", "Coasting", "Eclipsed"], 
        #         position = :rt,
        #         height = 0.8ppi,
        #         framevisible = false,
        #     )
        # end

        save(joinpath(dir_path, "transfer.pdf"), fig)
    end
    return nothing
end

function main()
    # Results files
    files = Dict(
        "gto_min_fuel" => joinpath(@__DIR__, "..", "data", "TAES", "GEO_MinFuel.jld2"),
        "gto_min_time" => joinpath(@__DIR__, "..", "data", "TAES", "GEO_MinTime.jld2"),
        "molniya_like_min_fuel" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLIke_MinFuel.jld2"),
        "molniya_like_min_time" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLIke_MinTime.jld2"),
        "molniya_like_np_min_fuel" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLikeNP_MinFuel.jld2"),
        "molniya_like_np_min_time" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLikeNP_MinTime.jld2"),
    )

    # Results to include in paper
    publication_cases = ("gto_min_fuel", "gto_min_time", "molniya_like_min_fuel", "molniya_like_min_time")

    for case in publication_cases
        # Load file
        qLawPs = load(files[case], "params")::QLawIndirectOptimization.qLawParams
        cache  = load(files[case], "cache")::QLawIndirectOptimization.QLawTransferCache

        # Call correct figure generation function
        dir = joinpath(@__DIR__, "..", "data", "TAES", "figures", "publication", case)
        if contains(case, "gto")
            #generate_gto_to_geo_figures(dir, cache, qLawPs)
        else
            generate_molniya_like_figures(dir, cache, qLawPs)
        end
    end
    return nothing
end

main()