using BayesHistogram, Test, Random

@testset "basic test BayesHistogram" begin
    x = randn(Xoshiro(1337), 5000)

    ref = [
        -3.704058531618975,
        -2.778065586707626,
        -2.133218533025471,
        -1.6480386516207801,
        -1.3140645798864168,
        -0.8575881683262867,
        -0.644610676860989,
        0.463534347245104,
        0.9542073946864447,
        1.2602857648643435,
        1.6456367564614374,
        1.9376034682498504,
        2.4043125317739733,
        2.932943071805272,
        3.590507153644763,
    ]
    @test all(bayesian_blocks(x, resolution = 100.0, min_counts = 2).edges .≈ ref)

    ref = [
        -3.704058531618975,
        -2.9198150710603223,
        -2.133218533025471,
        -1.3689581844499061,
        -0.6343551258945501,
        0.4495602463983795,
        1.1802613850539956,
        1.9231398688395287,
        2.6556538491763684,
        3.590507153644763,
    ]
    @test all(bayesian_blocks(x, resolution = 10.0, min_counts = 2).edges .≈ ref)

    ref = [
        -3.704058531618975,
        -2.3560907704049026,
        -1.8580502723199308,
        -1.3689581844499061,
        -0.8575881683262867,
        0.7468356102958359,
        1.2602857648643435,
        1.8727595838347098,
        2.372046297337469,
        3.590507153644763,
    ]
    @test all(bayesian_blocks(x, resolution = 15.0).edges .≈ ref)

    w = rand(Xoshiro(1337), 5000) .* 10 .+ 10
    ref = [
        -3.704058531618975,
        -2.133218533025471,
        -0.6740516063769897,
        1.747680727789755,
        3.590507153644763,
    ]
    @test all(bayesian_blocks(x, weights = w, resolution = 5.0).edges .≈ ref)
end



@testset "Jeffrey's Prior" begin
    max_N = 100
    Nv = 1:max_N
    prior = @. 1/sqrt(Nv)
    prior .= log.(prior./sum(prior))
    a_prior = Jeffrey(1.0)
    deltas = abs.(prior .- a_prior.(Nv, max_N))
    @test all(deltas .<= 1e-10)
end