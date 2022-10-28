using BayesHistogram, Test, Random, StableRNGs
const x = randn(StableRNG(1337), 5000)
const default_prior = Pearson(0.05)
naive_mean(f, x) = sum(f,x)/length(x)
naive_mean(x) = naive_mean(identity,x)

@testset "basic test BayesHistogram 1" begin
    ref = [
        -3.1578070050937224,
        -1.670203374306697,
        -1.2773998757939253,
        -0.9991635954562337,
        -0.8292472561339733,
        -0.6729400914817529,
        -0.5251263054226805,
        -0.3900877718108131,
        -0.23294418890701424,
        -0.1170237153552331,
        0.025508661047609786,
        0.14814423689292416,
        0.27664709881232963,
        0.403665425797619,
        0.5548106908002516,
        0.6970720044934977,
        0.8483378806816599,
        1.0630104033383305,
        1.3178383403513823,
        1.7117859538444158,
        3.8845391427592286,
    ]
    bl = bayesian_blocks(x, prior = default_prior, resolution = 100.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 2" begin

    ref = [
        -3.1578070050937224,
        -1.4787425372041283,
        -0.772606269519347,
        -0.06829690426062494,
        0.6361450014395109,
        1.3407574143055023,
        2.0454077110611504,
        3.8845391427592286,
    ]

    bl = bayesian_blocks(x, prior = default_prior, resolution = 10.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 3" begin

    ref = [
        -3.1578070050937224,
        -1.7812109359918114,
        -1.3028663292438303,
        -0.8327030469135661,
        -0.3629348545557658,
        0.10695053261983944,
        0.576474569553808,
        1.0460945821141636,
        1.5166659346731304,
        2.017887657637233,
        3.8845391427592286,
    ]
    bl = bayesian_blocks(x, prior = default_prior, resolution = 15.0).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 4" begin
    w = rand(StableRNG(1337), 5000) .* 10 .+ 10

    ref = [
        -3.1578070050937224,
        -1.4311770428250936,
        -0.02182252865867781,
        1.3869351312751794,
        3.8845391427592286,
    ]
    bl = bayesian_blocks(x, prior = default_prior, weights = w, resolution = 5.0).edges

    @test all(bl .≈ ref)

    ref = [
        -3.1578070050937224,
        -1.719768706521348,
        -1.3396850608727018,
        -1.0878481005264031,
        -0.878168482167752,
        -0.7301304287918617,
        -0.5718514075459757,
        -0.43089917727533195,
        -0.2845191851540201,
        -0.14365657389301362,
        -0.0026281281393301394,
        0.1386074703758985,
        0.2795515007275242,
        0.42044160300668515,
        0.5645770713245654,
        0.7055517425556006,
        0.8475238519063897,
        1.0630104033383305,
        1.3178383403513823,
        1.7117859538444158,
        3.8845391427592286,
    ]

    bl = bayesian_blocks(x, prior = default_prior, weights = w, resolution = 50.0).edges

    @test all(bl .≈ ref)
end

@testset "No Prior" begin
    a_prior = NoPrior()
    @test a_prior(0, 0, 0) == 0
end

@testset "Robustness" begin
    for pr in [Pearson(0.05), Geometric(0.05), Scargle(0.05), NoPrior(), BIC(), AIC(), HQIC()]
        tries = 100
        fails = 0
        okays = 0
        for i = 1:tries
            data = randn(StableRNG(1337 + 137 * i), 500)
            try
                bl = bayesian_blocks(data, prior = pr)
                okays += 1
            catch all
                fails += 1
            end
        end
        @test fails == 0
        @test okays == tries
    end
end

@testset "Integer inputs" begin
    input_int = round.(Int, randn(1000) * 100)
    input_float = float.(input_int)
    @test bayesian_blocks(input_int).counts == bayesian_blocks(input_float).counts
end

@testset "sanitize" begin
    xc = @. round(100 * x) / 100
    w = fill(1, 5000)
    w2 = w.^2
    idx = shuffle(StableRNG(1337), 1:5000)[1:4000]
    w[idx] .= 0
    w2[idx] .= 0
    xt, wt, wt2 = BayesHistogram.sanitize(xc, w, w2)
    @test all(xt .== unique(xt))
    @test all(wt .> 0)
    @test issorted(xt)
    @test sum(wt) == sum(w)
    @test sum(wt2) == sum(w2)
end

@testset "trivial datasets" begin
    @test bayesian_blocks([1.0]).counts[1] == 1
    @test bayesian_blocks(Float64[]).counts[1] == 0
end

@testset "error statistics" begin
    wh_dists = [
        rand(StableRNG(1338), length(x)), 
        one.(x),
        randexp(StableRNG(1338), length(x)),
    ]
    for wh in wh_dists
        wh2 = wh.^2
        bl = bayesian_blocks(x, weights = wh, prior = default_prior)
        m0 = 0
        m1 = zero.(bl.counts)
        m2 = zero.(bl.counts)
        for i in 1:500
            bx = rand(StableRNG(1337+i*137), x, length(x))
            sx,sw,sw2 = BayesHistogram.sanitize(bx,wh,wh2)
            cn = BayesHistogram.count_between_edges(bl.edges,sw,sx,false)
            m0+=1
            m1.+=cn
            m2.+=cn.^2
        end
        m1./=m0
        m2./=m0
        bs_err = sqrt.(m2 .- m1.^2)
        ratio_diff = abs.(bs_err ./ bl.error_counts .- 1)
        @test sort(ratio_diff)[end÷2] < 0.06
    end
end

@testset "error with custom sumw2" begin
    wh = rand(StableRNG(1338), length(x))
    wh_rebin2 = map(1:2:length(x)-1) do i
        sum(wh[i] + wh[i+1])
    end
    sumw2_rebin2 = map(1:2:length(x)-1) do i
        sum(wh[i]^2 + wh[i+1]^2)
    end
    x_rebin2 = x[begin:2:end]
    bl = bayesian_blocks(x, weights = wh, prior=Pearson(0.5))
    blwrong = bayesian_blocks(x_rebin2, weights = wh_rebin2, prior=Pearson(0.5))
    blright = bayesian_blocks(x_rebin2, weights = wh_rebin2, prior=Pearson(0.5),
                        sumw2 = sumw2_rebin2)
    @test blwrong.error_counts ≉ blright.error_counts atol=3
    @test bl.error_counts ≈ blright.error_counts atol=0.5
end

@testset "test severely rounded data" begin
    rng = StableRNG(123514)
    N1 = 3000
    N2 = 1000

    data = shuffle(
        rng,
        [
            randexp(rng, N1)
            randn(rng, N2) .* 0.02 .+ sqrt(1)
            randn(rng, N2) .* 0.02 .+ sqrt(2)
            randn(rng, N2) .* 0.04 .+ sqrt(4)
            randn(rng, N2) .* 0.08 .+ sqrt(8)
        ],
    )

    bl_noround = bayesian_blocks(data, prior = Pearson(0.2))
    @. data = round(data * 100) / 100
    bl_round = bayesian_blocks(data, prior = Pearson(0.2))
    @test sum(abs2,(bl_noround.edges .- bl_round.edges)) < 0.0005
    @test sum(abs2,(bl_noround.error_counts .- bl_round.error_counts)) < 0.2
end

@testset "to_pdf" begin
    bl = bayesian_blocks(x, prior=BIC())
    res = to_pdf(bl)
    @test res[1] == [-3.1578070050937224, -3.1578070050937224, -2.494702398745103, -2.494702398745103, -1.8298532494637834, -1.8298532494637834, -1.4561209721822812, -1.4561209721822812, -0.9536217995699808, -0.9536217995699808, 0.8030403893027767, 0.8030403893027767, 1.1970293913971268, 1.1970293913971268, 1.6398624165469307, 1.6398624165469307, 2.1005843440399232, 2.1005843440399232, 2.498821309397118, 2.498821309397118, 2.92140180959945, 2.92140180959945, 3.8845391427592286, 3.8845391427592286]
    @test res[2] == [0.0004845276479270787, 0.008143511518846264, 0.008143511518846264, 0.03399267341235206, 0.03399267341235206, 0.10970419867994981, 0.10970419867994981, 0.18666697402180743, 0.18666697402180743, 0.3594316562395907, 0.3594316562395907, 0.23909296832971377, 0.23909296832971377, 0.14316908721645752, 0.14316908721645752, 0.06511519901656145, 0.06511519901656145, 0.032141667181797555, 0.032141667181797555, 0.009465652101989552, 0.009465652101989552, 0.001453582943781236, 0.001453582943781236, 0.0004845276479270787]
end

@testset "estimate" begin
    bl = bayesian_blocks(x, prior = AIC())
    for fn in [one, identity, abs, abs2, cos, tanh]
        mu_b, sd_b = estimate(fn, bl)
        mu = naive_mean(fn,x)
        sd = sqrt(naive_mean(v->abs2(fn(v)-mu),x)/length(x))
        @test mu ≈ mu_b atol = max(sd,sd_b)
        @test sd ≈ sd_b atol = 0.02
    end
end