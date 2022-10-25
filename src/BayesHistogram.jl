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

function (w::NoPrior)(max_blocks, cnt_total, cnt_single)
    zero(cnt_single)
end

struct Geometric{T<:Real}
    gamma::T
end
function (w::Geometric)(max_blocks, cnt_total, cnt_single)
    # the normalisation constant can be omitted
    return -log(w.gamma)
end



struct Pearson{T<:Real}
    p::T
end
function (w::Pearson)(max_blocks, cnt_total, cnt_single)
    #                   unused
    return -cnt_total*(cnt_single/cnt_total - w.p)^2/w.p
end


struct Scargle{T<:Real}
    p0::T
end
function (w::Scargle)(max_blocks, cnt_total, cnt_single)
    #                              unused      unused
    C0 = 73.53
    C1 = -0.478
    log(C0 * w.p0 * max_blocks^C1) - 4.0
end

function bayesian_blocks(
    t::AbstractVector{T};
    weights::AbstractVector{W} = T[],
    prior = Pearson(0.05),
    resolution = Inf,
    min_counts::Integer = 0,
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
    change_points = Int[]
    sizehint!(change_points, N)
    ind = N + 1
    while true
        push!(change_points, ind)
        ind <= 1 && break
        ind = lasts[ind-1]
    end
    reverse!(change_points)
    edges = edges[change_points]

    # Evaluate densities and heights
    centers = @views(edges[begin:end-1] .+ edges[begin+1:end]) ./ 2
    counts = count_between_edges(edges, weights, t)
    total = sum(counts)
    widths = diff(edges)
    heights = counts ./ (total .* widths)
    return (; edges, counts, centers, widths, heights)
end

export bayesian_blocks, Pearson, Geometric, Scargle, NoPrior
end
