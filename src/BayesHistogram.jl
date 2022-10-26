"""BayesHistogram.jl
    main procedure: 
    bayesian_blocks(
        t::AbstractVector{T};
        weights::AbstractVector{W} = T[],
        prior = Pearson(0.05),
        resolution = Inf,
        min_counts::Integer = 0,
    )
"""
module BayesHistogram

function sort_sane(raw_x, raw_w)
    p = sortperm(raw_x)
    tmp_x = raw_x[p]
    tmp_w = raw_w[p]
    f = tmp_w .> 0
    tmp_x[f], tmp_w[f]
end

function sanitize(raw_x, raw_w)
    # todo: make more efficient
    if length(raw_w) == 0
        raw_w = one.(raw_x)
    end
    # first round:
    # - sort by x
    # - skip values with weights < 0
    x, w = sort_sane(raw_x, raw_w)
    # second round:
    # - merge entries with same x value
    m = fill(true, length(x))
    for i = 2:length(x)
        if x[i] == x[i-1] && m[i-1] && m[i] # merge
            m[i-1] = false
            w[i] += w[i-1]
            w[i-1] = 0
        end
    end
    return x[m], w[m]
end

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
    Geometric(x::T) where {T} =
        0 < x < 1 ? new{T}(x) : error("Gamma parameter must be between 0 and 1.")
end

function (w::Geometric)(max_blocks, cnt_total, cnt_single)
    # the normalisation constant can be omitted
    return -log(w.gamma / (1 - w.gamma))
end



struct Pearson{T<:Real}
    p::T
    Pearson(x::T) where {T} =
        0 < x < 1 ? new{T}(x) : error("probability parameter must be between 0 and 1.")
end
function (w::Pearson)(max_blocks, cnt_total, cnt_single)
    #                   unused
    return -cnt_total * (cnt_single / cnt_total - w.p)^2 / w.p
end


struct Scargle{T<:Real}
    p0::T
    Scargle(x::T) where {T} =
        0 < x < 1 ? new{T}(x) :
        error("false positive rate parameter must be between 0 and 1.")
end

function (w::Scargle)(max_blocks, cnt_total, cnt_single)
    #                              unused      unused
    C0 = 73.53
    C1 = -0.478
    log(C0 * w.p0 * max_blocks^C1) - 4.0
end

function build_blocks(t, edges, weights)
    centers = @views(edges[begin:end-1] .+ edges[begin+1:end]) ./ 2
    counts = count_between_edges(edges, weights, t)
    counts0 = count_between_edges(edges, one.(weights), t)
    counts2 = count_between_edges(edges, weights .^ 2, t)
    error_counts = sqrt.(counts2 .- (counts .^ 2) ./ counts0)
    total = sum(counts)
    widths = diff(edges)
    heights = counts ./ (total .* widths)
    error_heights = error_counts ./ (total .* widths)
    return (; edges, counts, centers, widths, heights, error_counts, error_heights)
end

function bayesian_blocks(
    t::AbstractVector{T};
    weights::AbstractVector{W} = T[],
    prior = Pearson(0.05),
    resolution = Inf,
    min_counts::Integer = 0,
) where {T<:Real,W<:Real}
    # copy and sort the arrays
    t, weights = sanitize(t, weights)
    
    # check trivial cases
    N = length(t)
    if N == 0
        return build_blocks(T[], T[-Inf, Inf], T[])
    elseif N == 1
        return build_blocks(t, T[t[1], t[1]], weights)
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
    return build_blocks(t, edges, weights)
end

export bayesian_blocks, Pearson, Geometric, Scargle, NoPrior
end
