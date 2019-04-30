using Test

@testset "Binning" begin
    include("BinningTests.jl")
    include("ProbabilityTests.jl")
    include("TagComparisonTests.jl")
end
