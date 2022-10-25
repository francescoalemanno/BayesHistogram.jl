using BayesHistogram, Test, Random
const x = randn(Xoshiro(1337), 5000)

@testset "basic test BayesHistogram 1" begin
    ref = [-3.704058531618975, -2.778065586707626, -2.133218533025471, -1.6480386516207801, -1.3140645798864168, -0.8575881683262867, 0.9542073946864447, 1.2602857648643435, 1.6456367564614374, 1.9376034682498504, 2.4043125317739733, 2.932943071805272, 3.590507153644763]

    bl = bayesian_blocks(x, resolution = 100.0, min_counts = 2).edges
    @show bl
    @test all(bl .≈ ref)
end
@testset "basic test BayesHistogram 2" begin

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
    bl = bayesian_blocks(x, resolution = 10.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end
@testset "basic test BayesHistogram 3" begin

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
    bl = bayesian_blocks(x, resolution = 15.0).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 4" begin
    w = rand(Xoshiro(1337), 5000) .* 10 .+ 10

    ref = [
        -3.704058531618975,
        -2.133218533025471,
        -0.6740516063769897,
        1.747680727789755,
        3.590507153644763,
    ]
    bl = bayesian_blocks(x, weights = w, resolution = 5.0).edges
    @test all(bl .≈ ref)
    
    ref = [-3.704058531618975, -2.778065586707626, -2.4631424510792215, -2.314878357736271, -2.1402068536709837, -1.9939563535400828, -1.8463434351082413, -1.6280845752313406, -1.3140645798864168, -1.089256530522866, -0.8575881683262867, -0.644610676860989, -0.49529999221146453, -0.3473371123214944, -0.06647876496714035, 0.09072830817794457, 0.29740423003699257, 0.463534347245104, 0.6864434760871314, 0.9542073946864447, 1.2602857648643435, 1.4805063922752275, 1.6456367564614374, 1.8727595838347098, 2.023652674002043, 2.289175115404521, 2.4649006788639394, 2.8689255970065712, 3.590507153644763]
    
    bl = bayesian_blocks(x, weights = w, resolution = 50.0).edges
    @test all(bl .≈ ref)
end

@testset "Jeffrey's Prior" begin
    max_N = 100
    Nv = 1:max_N
    prior = @. 1/sqrt(Nv)
    prior .= log.(prior./sum(prior))
    a_prior = Jeffrey(0.0)
    deltas = abs.(prior .- a_prior.(0, max_N, Nv))
    @test all(deltas .<= 1e-10)
end

@testset "No Prior" begin
    a_prior = NoPrior()
    @test a_prior(0,0,0) == 0
end

@testset "Robustness" begin
    tries = 100
    fails = 0
    okays = 0
    for i in 1:tries
        data = randn(Xoshiro(1337 + 137*i), 500)
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

