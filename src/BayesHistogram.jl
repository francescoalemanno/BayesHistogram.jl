"""BayesHistogram.jl
    main procedure: 
    function bayesian_blocks(
        t::AbstractVector{T};
        weights::AbstractVector{W} = one.(t),
        sumw2::AbstractVector{W} = abs2.(weights),
        prior = BIC(),
        resolution = Inf,
        min_counts::Real = 0,
    ) where {T<:Real,W<:Real}
"""
module BayesHistogram

include("priors.jl")

function sort_sane(raw_x, raw_w, raw_w2)
    p = sortperm(raw_x)
    tmp_x = raw_x[p]
    tmp_w = raw_w[p]
    tmp_w2 = raw_w2[p]
    f = tmp_w .> 0
    tmp_x[f], tmp_w[f], tmp_w2[f]
end

function sanitize(raw_x, raw_w, raw_w2)
    # first round:
    # - sort by x
    # - skip values with weights < 0
    x, w, w2 = sort_sane(raw_x, raw_w, raw_w2)
    # second round:
    # - merge entries with same x value
    m = fill(true, length(x))
    for i = 2:length(x)
        if x[i] == x[i-1] && m[i-1] && m[i] # merge
            m[i-1] = false
            w[i] += w[i-1]
            w2[i] += w2[i-1]
            w[i-1] = 0
            w2[i-1] = 0
        end
    end
    return x[m], w[m], w2[m]
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

function build_blocks(t, edges, weights, weights2)
    centers = @views(edges[begin:end-1] .+ edges[begin+1:end]) ./ 2
    counts = count_between_edges(edges, weights, t)
    counts2 = count_between_edges(edges, weights2, t)
    total = sum(counts)
    error_counts = sqrt.(max.(counts2 .- counts.^2 ./ total, 0))
    widths = diff(edges)
    heights = counts ./ (total .* widths)
    error_heights = error_counts ./ (total .* widths)
    return (; edges, counts, centers, widths, heights, error_counts, error_heights)
end

"""
    function bayesian_blocks(
        datas::AbstractVector{T};
        weights::AbstractVector{W} = one.(t),
        sumw2::AbstractVector{W} = abs2.(weights),
        prior = BIC(),
        resolution = Inf,
        min_counts::Real = 0,
    ) where {T<:Real,W<:Real}

- `datas`: Observations
- `weights`: sample weight of each observation
- `sumw2`: sum of weight^2 in each observation, this is particularly useful
when this algorihtm is used for re-binning of already made histograms where
the `sumw2` for each bin is different from `weight^2` of each bin.
- `prior`: choose from `NoPrior`, `Pearson`, `Geometric`, `BIC`, `AIC`, `HQIC`
- resolution: handles on how fine we count along the `datas axis
- min_counts: minimum sum of weights of a block that can be splitted.

"""
function bayesian_blocks(
    t::AbstractVector{T};
    weights::AbstractVector{W} = one.(t),
    sumw2::AbstractVector{W} = abs2.(weights),
    prior = BIC(),
    resolution = Inf,
    min_counts::Real = 0,
) where {T<:Real,W<:Real}
    # copy and sort the arrays
    t, weights, sumw2 = sanitize(t, weights, sumw2)
    
    # check trivial cases
    N = length(t)
    if N == 0
        return build_blocks(T[], T[-Inf, Inf], T[], sumw2)
    elseif N == 1
        return build_blocks(t, T[t[1], t[1]], weights, sumw2)
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
    return build_blocks(t, edges, weights, sumw2)
end

export bayesian_blocks, Pearson, Geometric, Scargle, NoPrior, BIC, AIC, HQIC
end
