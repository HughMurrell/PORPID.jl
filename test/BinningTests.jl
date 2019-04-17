include("../src/include_all.jl")
using Test

@testset "Binning" begin
  @testset "IO and Setup" begin
    include("TestConfig.jl")
    global cfg = example_config()
    include("TestSequences.jl")
  end

end
