#has been run only with julia 1.8+
using Plots, Random, Statistics
using BayesHistogram
gr()


function fake_data_pdf(v)
    # a fake statistical model with peaks and exponential background
    N1 = 3000
    tN = N1
    x = N1*exp(-v)
    for i in 1:5
        N = round(Int64,0.4*N1/sqrt(i))
        tN += N
        std = i*0.02
        mu = i*0.7
        x += N*exp(-abs2(v-mu)/(2*std^2))/sqrt(2*pi*std^2)
    end
    return x/tN
end

function fake_data()
    # sampler for the fake model above
    rng = Xoshiro(123513)
    N1 = 3000
    x = randexp(rng, N1)
    for i in 1:5
        N = round(Int64,0.4*N1/sqrt(i))
        std = i*0.02
        mu = i*0.7
        append!(x, randn(rng, N) .* std .+ mu)
    end
    @. x = round(x * 200) / 200
    return x
end

let plot_mode = true
    x = fake_data()
    stephist(x, normalize = :pdf, lw=1.7, label = "classical", legend=:topright)

    b1 = bayesian_blocks(x, prior=HQIC())
    if plot_mode
        # this uses BayesHistogram.jl "to_pdf"
        xvals, epdf = to_pdf(b1)
        plot!(xvals, epdf, lw=2, label = "bayes-hist")
    else
        # this uses Plots.jl "stephist"
        stephist!(x, bin = b1.edges, normalize = :pdf,lw=2, label = "bayes-hist")    
    end
    lx = 0:0.001:8
    pdf = fake_data_pdf.(lx)
    plot!(lx,pdf, lw = 1.5,color="black", linestyle = :dash, label = "density")
    xlims!(0,4.5)
    savefig("plot3.png")
    plot!()
end


let
    gr()
    x = fake_data()
    b = bayesian_blocks(x, prior=HQIC())
    P = []

    Hs = [
        ("Bayes Histogram (def. mode)", 0),
        ("Sturges", :sturges),
        ("Rice", :rice),
        ("Sqrt", :sqrt),
    ]
    i = 1
    for (lab, bin) in Hs        
        pl = plot(yaxis = :log,grid = false,)
        if i == 1 
            X,Y = to_pdf(b)
            plot!(pl, X, Y,label = :none)
        else
            stephist!(pl,
                x,
                label = :none,
                bins = bin,
                yaxis = :log,
                normalize = :pdf,
                lw = 1,
            )
        end
        title!(pl, lab, titlefont = font(12))
        push!(P, pl)
        i += 1
    end

    p = plot!(P..., layout = @layout([a b; c d]), size = (500, 300) .* 1.3)
    savefig(p, "plot.png")
    p
end


let
    gr()
    rng = Xoshiro(1231)
    x = collect(-4:0.01:4)
    w = round.(Int,exp.(.-x.^2).*30 .+ rand(rng, length(x))) .+ 0.5
    b = bayesian_blocks(x, weights = w)
    scatter(x, shuffle(rng,LinRange(0,0.5, length(x))), ms=sqrt.(w),alpha=0.5, label="noisy weighted data",markercolor="red",markerstrokewidth=0)
    dx, dy = to_pdf(b)

    plot!(
        dx, dy,
        color = "black",
        label = "bayes-hist",
        lw = 2.5,
    )

    scatter!(
        b.centers,
        b.heights,
        yerr = b.error_heights,
        label = :none,
        color = "black",
        lw = 2.5,
    )

    savefig("plot2.png")
    plot!()
end
