
using QLawIndirectOptimization
using JLD2, FileIO
using CairoMakie

function generate_gto_to_geo_figures(
    dir_path, cache::QLawTransferCache, ps::qLawParams;
)
    axes = (SA[1,2], SA[2,3])
    return nothing
end

function generate_molniya_like_figures(
    dir_path, cache::QLawTransferCache, ps::qLawParams;
)

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

    for case in eachindex(publication_cases)
        # Load file
        qLawPs = load(files[case], "params")::QLawIndirectOptimization.qLawParams
        cache  = load(files[case], "cache")::QLawIndirectOptimization.QLawTransferCache

        # Call correct figure generation function
        dir = joinpath(@__DIR__, "..", "data", "TAES", "figures", "publication", case)
        if contains(case, "gto")
            generate_gto_to_geo_figures(dir, cache, qLawPs)
        else
            generate_molniya_like_figures(dir, cache, qLawPs)
        end
    end
    return nothing
end

main()