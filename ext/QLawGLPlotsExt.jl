
module QLawGLPlotsExt

using QLawIndirectOptimization
import GLMakie as GM

function plot_transfer_gl(
    cache::QLawIndirectOptimization.QLawTransferCache, 
    ps::QLawIndirectOptimization.qLawParams,
)
    fig = GM.Figure();
    ax = GM.Axis3(fig[1, 1], aspect = :data)
    xs = zeros(length(cache.times))
    ys = zeros(length(cache.times))
    zs = zeros(length(cache.times))
    for i in eachindex(cache.times)
        cart = cache.states[i]
        xs[i] = cart[1]
        ys[i] = cart[2]
        zs[i] = cart[3]
    end
    GM.lines!(ax, xs, ys, zs)
    return fig
end

end