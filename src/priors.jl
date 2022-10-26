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

