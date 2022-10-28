
function to_pdf(h::BHist{T}; lb = minimum(h.heights)/3) where T
    X = T[]
    Y = T[]
    push!(X, h.edges[begin])
    push!(Y, lb)
    for i in eachindex(h.heights)
        push!(X, h.edges[i])
        push!(X, h.edges[i+1])
        push!(Y, h.heights[i])
        push!(Y, h.heights[i])
    end
    push!(X, h.edges[end])
    push!(Y, lb)
    return X, Y
end

function estimate(f::Function, h::BHist{T}) where T
    function element(i::Int)
        x = h.centers[i]
        ht = h.heights[i]
        wd = h.widths[i]
        A = ht * wd
        eA = h.error_heights[i] * wd
        fv = f(x)
        ii = fv*A
        eii = abs2(fv*eA)
        (ii, eii)::Tuple{T,T}
    end
    s = element(1)
    for i in 2:length(h.centers)
        s = s .+ element(i)
    end
    s[1], sqrt(s[2])
end

export to_pdf, estimate