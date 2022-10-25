using Plots, Random, Statistics
using BayesHistogram
let
    gr()
    rng = Xoshiro(123514)
    N1 = 3000
    N2 = 1000

    x = shuffle(rng,[
        randexp(rng,N1); 
        randn(rng,N2).*0.02 .+ sqrt(1); 
        randn(rng,N2).*0.02 .+ sqrt(2); 
        randn(rng,N2).*0.04 .+ sqrt(4); 
        randn(rng,N2).*0.08 .+ sqrt(8);
    ])

    b = bayesian_blocks(x)
    P = []

    edg_equi_area = quantile(x, range(0,1,length=ceil(Int, 2*length(x)^(1.88/5))))
    Hs = [
        ("Bayes Histogram", b.edges),
        ("EquiArea", edg_equi_area),
        ("Rice", :rice),
        ("Sqrt", :sqrt),
    ]

    for (lab, bin) in Hs
        pl = stephist(x, label=:none, bins=bin, yaxis=:log, normalize=:density, lw=1, grid=false)
        title!(pl, lab, titlefont = font(12))
        push!(P, pl)
    end

    p = plot!(P..., layout = @layout([a b;c d]), size=(500,300).*1.3)
    savefig(p, "plot.png")
    p
end


using Trapz
let
    gr()
    rng = Xoshiro(123514)
    x = collect(range(-4,4,length = 300))
    w = ceil.(rand(rng,length(x)).*100 .+ @. ( exp(-x^2/2)*500 .+ 1))
    #in this case let's use the non-informative prior
    b = bayesian_blocks(x, weights=w, prior = Jeffrey(3.0))
    pdf = w./trapz(x,w)
    plot(x,pdf,color="red", label="noisy weighted obs")
    stephist!(x,bins = b.edges, weights=w, normalize=true, color="black", yaxis=:log,label="bayeshist",lw=2)
    savefig("plot2.png")
    plot!()
end