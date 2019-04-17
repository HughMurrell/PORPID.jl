include("../src/include_all.jl")
using Test

@testset "Binning" begin
    include("IOTests.jl")
    #include("test_fasta.jl")
    #include("test_fastq.jl")
end
