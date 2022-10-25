using BayesHistogram, Test, Random, StableRNGs
const x = randn(StableRNG(1337), 5000)

@testset "basic test BayesHistogram 1" begin
    ref = [-3.1578070050937224, -1.670203374306697, -1.2773998757939253, -0.9991635954562337, -0.8292472561339733, -0.6729400914817529, -0.5251263054226805, -0.3900877718108131, -0.23294418890701424, -0.1170237153552331, 0.025508661047609786, 0.14814423689292416, 0.27664709881232963, 0.403665425797619, 0.5548106908002516, 0.6970720044934977, 0.8483378806816599, 1.0630104033383305, 1.3178383403513823, 1.7117859538444158, 3.8845391427592286]
    bl = bayesian_blocks(x, resolution = 100.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 2" begin

    ref = [-3.1578070050937224, -1.4787425372041283, -0.772606269519347, -0.06829690426062494, 0.6361450014395109, 1.3407574143055023, 2.0454077110611504, 3.8845391427592286]

    bl = bayesian_blocks(x, resolution = 10.0, min_counts = 2).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 3" begin

    ref =  [-3.1578070050937224, -1.7812109359918114, -1.3028663292438303, -0.8327030469135661, -0.3629348545557658, 0.10695053261983944, 0.576474569553808, 1.0460945821141636, 1.5166659346731304, 2.017887657637233, 3.8845391427592286]
    bl = bayesian_blocks(x, resolution = 15.0).edges
    @test all(bl .≈ ref)
end

@testset "basic test BayesHistogram 4" begin
    w = rand(StableRNG(1337), 5000) .* 10 .+ 10

    ref = [-3.1578070050937224, -1.4311770428250936, -0.02182252865867781, 1.3869351312751794, 3.8845391427592286]
    bl = bayesian_blocks(x, weights = w, resolution = 5.0).edges
    
    @test all(bl .≈ ref)

    ref = [-3.1578070050937224, -1.719768706521348, -1.3396850608727018, -1.0878481005264031, -0.878168482167752, -0.7301304287918617, -0.5718514075459757, -0.43089917727533195, -0.2845191851540201, -0.14365657389301362, -0.0026281281393301394, 0.1386074703758985, 0.2795515007275242, 0.42044160300668515, 0.5645770713245654, 0.7055517425556006, 0.8475238519063897, 1.0630104033383305, 1.3178383403513823, 1.7117859538444158, 3.8845391427592286]

    bl = bayesian_blocks(x, weights = w, resolution = 50.0).edges

    @test all(bl .≈ ref)
end

@testset "No Prior" begin
    a_prior = NoPrior()
    @test a_prior(0, 0, 0) == 0
end

@testset "Robustness" begin
    for pr in [Pearson(0.05), Geometric(0.05), Scargle(0.05), NoPrior()]
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
    @test bayesian_blocks(input_int) == bayesian_blocks(input_float)
end
