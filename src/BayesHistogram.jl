"""BayesHistogram.jl
    main procedure: 
    bayesian_blocks(
        t::AbstractVector{T};
        weights::AbstractVector{W}=T[],
        prior = Scargle(T(0.05)),
        resolution = T(Inf),
        min_counts::Integer = -1
    )
"""
module BayesHistogram
function count_between_edges(edges, weights, observations, shift::Bool = false)
    i = 1
    out = zeros(length(edges) - 1 + shift)
    for (el, w) in zip(observations, weights)
        while !(edges[i] <= el <= edges[i+1])
            i += 1
        end
        out[i+shift] += w
    end
    return out
end

struct NoPrior end

function (w::NoPrior)(cnt_blocks, cnt_total, cnt_single)
    zero(cnt_single)
end

struct Jeffrey{T<:Real}
    penalty::T
end

function (w::Jeffrey)(cnt_blocks, cnt_total, cnt_single)
    #                   unused
    C0 = -0.020833333333333332
    C1 = -0.730177254404794
    max(1 + w.penalty, zero(w.penalty)) * log(
        (cnt_total * sqrt(cnt_total) / 2) /
        (sqrt(cnt_single) * (C0 + cnt_total * (1 / 4 + C1 * sqrt(cnt_total) + cnt_total))),
    )
end

struct Scargle{T<:Real}
    p0::T
end
function (w::Scargle)(cnt_blocks, cnt_total, cnt_single)
    #                              unused      unused
    C0 = 73.53
    C1 = -0.478
    log(C0 * w.p0 * cnt_blocks^C1) - 4.0
end

function bayesian_blocks(
    t::AbstractVector{T};
    weights::AbstractVector{W} = T[],
    prior = Scargle(0.05),
    resolution = Inf,
    min_counts::Integer = -1,
) where {T<:Real,W<:Real}
    N = length(t)
    # copy and sort the arrays
    perm = sortperm(t)
    t = t[perm]

    if length(weights) == 0
        weights = ones(T, N)
    else
        weights = weights[perm]
    end

    if min_counts <= -1
        min_counts = ceil(Int64, sqrt(sum(weights)) / 2)
    end
    # create cell edges
    edges = [t[begin]; @views(t[begin+1:end] .+ t[begin:end-1]) ./ 2; t[end]]

    # make cumulative weights
    wh_in_edge = count_between_edges(edges, weights, t, true)
    wh_in_edge .= cumsum(wh_in_edge)

    #=
    # this was for testing the weights distribution
    for i in 1:N, j in i:N
        c1 = sum(@views(weights[i:j])) == wh_in_edge[j+1]-wh_in_edge[i]
        c2 = edges[i] <= t[i] <= edges[j+1]
        c3 = edges[i] <= t[j] <= edges[j+1]
        if !(c1 || c2 || c3)
            error("wrong $i $j $N, $c1 $c2 $c3")
        end
    end
    =#

    # arrays needed for the iteration
    best = zeros(N)
    lasts = zeros(Int, N)

    extent = t[end] - t[begin]
    # by default dt is 0 because resolution is Inf
    dt = max(abs(extent / resolution), zero(T))
    # Start with first data cell; add one cell at each iteration
    @inbounds for Q = 1:N
        fit_max = -Inf
        i_max = 0
        for i = 1:Q
            cnt_in_range = wh_in_edge[Q+1] - wh_in_edge[i]
            cnt_in_range < min_counts && break
            width = edges[Q+1] - edges[i]
            width <= dt && break

            fitness =
                cnt_in_range * (log(cnt_in_range / width)) +
                prior(N, wh_in_edge[end], cnt_in_range)
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
    centers = @views(edges[begin:end-1] .+ edges[begin+1:end]) ./ 2
    counts = count_between_edges(edges, weights, t)
    total = sum(counts)
    widths = diff(edges)
    heights = counts ./ (total .* widths)
    return (; edges, counts, centers, widths, heights)
end

export bayesian_blocks, Jeffrey, Scargle, NoPrior
end
