include("../src/include_all.jl")
using Test

@testset "Binning" begin
  @testset "IO and Setup" begin
    include("TestConfig.jl")
    global cfg = example_config()
    include("TestSequences.jl")
  end

  function checkBins(source_file_name, binned_template, binned_tag, output_sequence, score)
      sequence_name, authority_template, authority_tag = split(FASTQ.identifier(output_sequence), "-")
      @test authority_template == binned_template.name
      @test authority_tag == binned_tag
  end

  @testset "Processing basic fastq" begin
    process_file("test_data/basic.fastq", cfg, checkBins)
  end
  @testset "Processing more than basic fastq" begin
    process_file("test_data/more_than_basic.fastq", cfg, checkBins)
  end
end
