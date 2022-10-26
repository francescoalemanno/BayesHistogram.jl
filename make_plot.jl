#has been run only with julia 1.8+
using Plots, Random, Statistics
using BayesHistogram
let
    gr()
    rng = Xoshiro(123514)
    N1 = 3000
    N2 = 1000

    x = shuffle(
        rng,
        [
            randexp(rng, N1)
            randn(rng, N2) .* 0.02 .+ sqrt(1)
            randn(rng, N2) .* 0.02 .+ sqrt(2)
            randn(rng, N2) .* 0.04 .+ sqrt(4)
            randn(rng, N2) .* 0.08 .+ sqrt(8)
        ],
    )
    @. x = round(x * 100) / 100
    b = bayesian_blocks(x, prior = BIC())
    P = []

    edg_equi_area = quantile(x, range(0, 1, length = ceil(Int, 2 * length(x)^(1.88 / 5))))
    Hs = [
        ("Bayes Histogram", b.edges),
        ("EquiArea", edg_equi_area),
        ("Rice", :rice),
        ("Sqrt", :sqrt),
    ]
    i = 1
    for (lab, bin) in Hs
        pl = stephist(
            x,
            label = :none,
            bins = bin,
            yaxis = :log,
            normalize = :pdf,
            lw = 1,
            grid = false,
        )
        i == 1 && scatter!(pl,
            b.centers,
            b.heights,
            yerr = b.error_heights,
            label = :none,
            lw = 2.0,
            markersize= 2.0
        )
        title!(pl, lab, titlefont = font(12))
        push!(P, pl)
        i += 1
    end

    p = plot!(P..., layout = @layout([a b; c d]), size = (500, 300) .* 1.3)
    savefig(p, "plot.png")
    p
end


using Trapz
let
    gr()
    rng = Xoshiro(12351)
    x = collect(range(-4, 4, length = 300))
    w = ceil.((rand(rng, length(x)) .* 3 .+ @. (exp(-x^2 / 2) * 15)))
    pdf = w ./ trapz(x, w)

    b = bayesian_blocks(x, weights = w, prior = AIC())

    plot(x, pdf, color = "red", label = "noisy weighted obs", lw = 0.5)
    stephist!(
        x,
        bins = b.edges,
        weights = w,
        normalize = :pdf,
        color = "black",
        label = "bayeshist",
        lw = 2,
        yaxis=:log,
    )
    scatter!(
        b.centers,
        b.heights,
        yerr = b.error_heights,
        label = :none,
        color = "black",
        lw = 2,
    )

    savefig("plot2.png")
    plot!()
end
