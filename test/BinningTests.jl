using Test

@testset "Binning" begin
  @testset "IO and Setup" begin
    include("TestConfig.jl")
    global cfg = example_config()
    include("TestSequences.jl")
  end

  function check_bins(source_file_name, binned_template, binned_tag, output_sequence, score)
      sequence_name, authority_template, authority_tag = split(FASTQ.identifier(output_sequence), "-")
      @test authority_template == binned_template.name
      @test authority_tag == binned_tag
  end

  @testset "Processing basic fastq" begin
    extract_tags_from_file("test_data/basic.fastq", cfg, check_bins)
  end
  @testset "Processing more than basic fastq" begin
    extract_tags_from_file("test_data/more_than_basic.fastq", cfg, check_bins)
  end
end
