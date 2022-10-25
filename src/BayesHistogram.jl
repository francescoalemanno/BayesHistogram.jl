"""BayesHistogram.jl
    main procedure: 
    bayesian_blocks(
        t::AbstractVector{T};
        p0::T = T(5) / T(100),
        resolution::T = T(Inf),
        min_counts::Int = ceil(Int64, sqrt(length(t))/2)
    )
"""
module BayesHistogram
function count_between_edges(edges,weights,observations, shift::Bool = false)
    i = 1
    out = zeros(length(edges) - 1 + shift)
    for (el,w) in zip(observations,weights)
        while !(edges[i] <= el < edges[i+1] || el == edges[end])
            i += 1
        end
        out[i+shift] += w
    end
    return out
end
function bayesian_blocks(
    t::AbstractVector{T};
    weights::AbstractVector{W}=T[],
    p0 = T(5) / T(100),
    resolution = oftype(p0, Inf),
    min_counts::Integer = -1
) where {T<:Real, W<:Real}
    N = length(t)
    # copy and sort the arrays
    perm = sortperm(t)
    t = t[perm]
    
    if length(weights) == 0
        weights = ones(T,N)
    else
        weights = weights[perm]
    end

    if min_counts <= -1
        min_counts = ceil(Int64, sqrt(sum(weights))/2)
    end
    # create cell edges
    edges = [t[begin]; @views(t[begin+1:end] .+ t[begin:end-1]) ./ 2; t[end]]

    # make cumulative weights
    wh_in_edge = count_between_edges(edges, weights, t, true)
    wh_in_edge .= reverse(cumsum(reverse(wh_in_edge)))
    # arrays needed for the iteration
    best = zeros(N)
    lasts = zeros(Int, N)

    extent = t[end] - t[begin]
    # by default dt is 0 because resolution is Inf
    dt = max(abs(extent / resolution), zero(T))
    lp0 = log(p0) + 0.2976934862081313

    # Start with first data cell; add one cell at each iteration
    @inbounds for Q = 1:N
        fit_max = -Inf
        i_max = 0
        for i = 1 : Q
            cnt_in_range = wh_in_edge[i] - wh_in_edge[Q+1]
            cnt_in_range < min_counts && break
            width = edges[Q+1] - edges[i]
            width <= dt && break

            fitness = cnt_in_range * (log(cnt_in_range / width)) + lp0 - 0.478 * log(cnt_in_range)
            if i > 1
                fitness += best[i-1]
            end
            if fitness > fit_max
                fit_max = fitness
                i_max = i
            end
        end

        # find the max of the fitness: this is the Q^th changepoint
        lasts[Q] = i_max
        best[Q] = fit_max
    end
    # Recover changepoints by iteratively peeling off the last block
    change_points = zeros(Int, N)
    i_cp = N + 1
    ind = N + 1
    while true
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 1
            break
        end
        ind = lasts[ind-1]
    end
    change_points = change_points[i_cp:end]
    edges = edges[change_points]
    
    # Evaluate densities and heights
    centers =  @views(edges[begin:end-1] .+ edges[begin+1:end]) ./ 2
    counts = count_between_edges(edges, weights, t)
    total = sum(counts)
    widths = diff(edges)
    heights = counts ./ (total .* widths)
    return (; edges, counts, centers, widths, heights)
end

export bayesian_blocks, count_between_edges
end
