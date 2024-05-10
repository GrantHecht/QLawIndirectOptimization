

using QLawIndirectOptimization
using JLD2, FileIO

function main(transfer_str::String)
    # Results files
    files = Dict(
        "gto_min_fuel" => joinpath(@__DIR__, "..", "data", "TAES", "GEO_MinFuel.jld2"),
        "gto_min_time" => joinpath(@__DIR__, "..", "data", "TAES", "GEO_MinTime.jld2"),
        "molniya_like_min_fuel" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLIke_MinFuel.jld2"),
        "molniya_like_min_time" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLIke_MinTime.jld2"),
        "molniya_like_np_min_fuel" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLikeNP_MinFuel.jld2"),
        "molniya_like_np_min_time" => joinpath(@__DIR__, "..", "data", "TAES", "MolniyaLikeNP_MinTime.jld2"),
    )

    # Load file
    qLawPs = load(files[transfer_str], "params")::QLawIndirectOptimization.qLawParams
    cache  = load(files[transfer_str], "cache")::QLawIndirectOptimization.QLawTransferCache

    # Precompuse results for output
    Δm = qLawPs.m0 - cache.states[end][7]
    nrev = cache.time_at_steps[end] / (2*pi)

    # Summarize results
    println("Transfer time:     $(cache.times[end]) days")
    println("Revolutions:       $nrev")
    println("Mass consumed:     $Δm kg")
    println("Weights:           [$(qLawPs.oeW[1]), $(qLawPs.oeW[2]), $(qLawPs.oeW[3]), $(qLawPs.oeW[4]), $(qLawPs.oeW[5])]")
    println("Eff. tolerance:    $(qLawPs.ηr)")

    # Generate figures
    generate_publication_figures(
        joinpath(@__DIR__, "..", "data", "TAES", "figures", "publication", transfer_str), cache, qLawPs
    )
end

gto_str_list = ("gto_min_time", "gto_min_fuel")
molniya_str_list = ("molniya_like_min_time", "molniya_like_np_min_time", "molniya_like_min_fuel", "molniya_like_np_min_fuel")
for str in molniya_str_list
    main(str)
end