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

function bayesian_blocks(
    t::AbstractVector{T};
    p0::T = T(5) / T(100),
    resolution::T = T(Inf),
    min_counts::Int = ceil(Int64, sqrt(length(t))/2)
) where {T<:AbstractFloat}
    # copy and sort the array
    t = sort(t)
    N = length(t)

    # create cell edges
    edges = [t[1]; (t[2:end] .+ t[1:end-1]) ./ 2; t[end]]
    block_length = t[end] .- edges

    # arrays needed for the iteration
    best = zeros(N)
    last = zeros(Int, N)

    extent = t[end] - t[1]
    dt = max(abs(extent / resolution), zero(extent))
    lp0 = log(p0) + 0.2976934862081313

    # Start with first data cell; add one cell at each iteration
    for Q = 1:N
        fit_max = -Inf
        i_max = 0
        for i = 1:Q
            width = block_length[i] - block_length[Q+1]
            width <= dt && continue
            cnt_in_range = Q + 1 - i
            cnt_in_range < min_counts && continue

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
        last[Q] = i_max
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
        ind = last[ind-1]
    end
    change_points = change_points[i_cp:end]
    edges = edges[change_points]

    # Evaluate densities and heights
    E = length(edges)
    counts = zeros(Int, E - 1)
    i = 1
    for el in t
        @label retry
        if edges[i] <= el < edges[i+1] || el == edges[end]
            counts[i] += 1
        else
            i += 1
            @goto retry
        end
    end
    total = sum(counts)
    heights = counts ./ (total .* diff(edges))
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    widths = diff(edges)
    return (; edges, counts, centers, widths, heights)
end

export bayesian_blocks
end