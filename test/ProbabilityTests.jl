using Test
using PORPID

@testset "Testing probabilities" begin
    @testset "Phred score to prob" begin
        # From the table on http://drive5.com/usearch/manual/quality_score.html
        @test phred_score_to_prob(0) ≈ 0.0
        @test phred_score_to_prob(10) ≈ 0.9
        @test phred_score_to_prob(20) ≈ 0.99
        @test isapprox(phred_score_to_prob(25), (1.0 - 0.00316); atol=1e-5)
        @test isapprox(phred_score_to_prob(42)  , (1.0 - 0.00006); atol=1e-5)
    end

    @testset "Probability of getting x given observation" begin
        @test prob(DNA_A, DNA_A, 0.9) ≈ 0.9
        @test prob(DNA_C, DNA_A, 0.9) ≈ 0.1/3
        @test prob(DNA_M, DNA_A, 1.0) ≈ 0.5
        #@test prob(DNA_M, DNA_A, 0.9) ≈ 0.9 + 0.1/3
        #@test prob(DNA_M, DNA_T, 0.9) ≈ 0.1*2/3
    end
end
