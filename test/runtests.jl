using BayesHistogram, Test, Random, StableRNGs
const x = randn(StableRNG(1337), 5000)

@testset "basic test BayesHistogram 1" begin
    ref = [
        -3.1578070050937224,
        -2.627589003399208,
        -2.1567159913141696,
        -1.809553033077446,
        -1.4561209721822812,
        -0.9972606510660637,
        -0.5753838312589055,
        0.8030403893027767,
        1.1970293913971268,
        1.6398624165469307,
        2.1005843440399232,
        2.498821309397118,
        2.92140180959945,
        3.8845391427592286,
    ]

    bl = bayesian_blocks(x, resolution = 100.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end
@testset "basic test BayesHistogram 2" begin

    ref = [
        -3.1578070050937224,
        -2.4424726736043123,
        -1.719768706521348,
        -0.9972606510660637,
        0.9449575195349826,
        1.6625393156205948,
        2.498821309397118,
        3.8845391427592286,
    ]

    bl = bayesian_blocks(x, resolution = 10.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end
@testset "basic test BayesHistogram 3" begin

    ref = [
        -3.1578070050937224,
        -2.1567159913141696,
        -1.670203374306697,
        -0.9972606510660637,
        -0.23294418890701424,
        0.3779185031870537,
        0.9593283778530831,
        1.4851079736094,
        1.965053158543435,
        2.4534081553862954,
        3.8845391427592286,
    ]
    bl = bayesian_blocks(x, resolution = 15.0).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 4" begin
    w = rand(StableRNG(1337), 5000) .* 10 .+ 10

    ref = [
        -3.1578070050937224,
        -1.5449036608365123,
        1.0655876974876393,
        2.474272965955527,
        3.8845391427592286,
    ]
    bl = bayesian_blocks(x, weights = w, resolution = 5.0).edges

    @test all(bl .≈ ref)

    ref = [
        -3.1578070050937224,
        -2.627589003399208,
        -2.4424726736043123,
        -2.1567159913141696,
        -1.836892002440313,
        -1.672177926449848,
        -1.4561209721822812,
        -1.2629080974042246,
        -0.9991635954562337,
        -0.7973904890758701,
        -0.5753838312589055,
        -0.4286389710947635,
        -0.23354006071343802,
        -0.08731284637243833,
        0.10695053261983944,
        0.3779185031870537,
        0.5255797020151315,
        0.8030403893027767,
        0.9449575195349826,
        1.1970293913971268,
        1.5203188316842984,
        1.7117859538444158,
        1.9626991729085734,
        2.104886260600142,
        2.2470590206020926,
        2.498821309397118,
        2.8530173295282495,
        3.8845391427592286,
    ]

    bl = bayesian_blocks(x, weights = w, resolution = 50.0).edges
    @show bl

    @test all(bl .≈ ref)
end

@testset "Jeffrey's Prior" begin
    max_N = 100
    Nv = 1:max_N
    prior = @. 1 / sqrt(Nv)
    prior .= log.(prior ./ sum(prior))
    a_prior = Jeffrey(0.0)
    deltas = abs.(prior .- a_prior.(0, max_N, Nv))
    @test all(deltas .<= 1e-10)
end

@testset "No Prior" begin
    a_prior = NoPrior()
    @test a_prior(0, 0, 0) == 0
end

@testset "Robustness" begin
    tries = 100
    fails = 0
    okays = 0
    for i = 1:tries
        data = randn(StableRNG(1337 + 137 * i), 500)
        try
            bl = bayesian_blocks(data)
            okays += 1
        catch all
            fails += 1
        end
    end
    @test fails == 0
    @test okays == tries
end


@testset "Integer inputs" begin
    input_int = round.(Int, randn(1000) * 100)
    input_float = float.(input_int)
    @test bayesian_blocks(input_int) == bayesian_blocks(input_float)
end
