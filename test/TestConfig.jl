using Test
using PORPID

Base.:(==)(x::Template, y::Template) = x.name == y.name && x.reference == y.reference

function example_config()
    templates = Dict()
    templates["GATTACA"] = "GATTACAGATTACAnnnnnnGATTACA*"
    templates["TAGACAT"] = "TAGACATTAGACATnnnnnnTAGACATATAT*"
    cfg = Configuration()
    cfg.files = ["test_data/basic.fastq", "test_data/more_than_basic.fastq"]
    cfg.filetype = fastq
    cfg.start_inclusive = 1
    #cfg.end_inclusive = ceil(Int, maximum([length(x) for x in values(templates)]) * 1.2)
    cfg.end_inclusive = 39
    cfg.max_allowed_errors = 3
    cfg.try_reverse_complement = true
    for (name, reference) in templates
        push!(cfg.templates, Template(name, reference))
    end
    return cfg
end

@testset "Creating config" begin

  # Test default config
  @testset "Default config doesn't have contradictory start and end points" begin
    default_cfg = Configuration()
    @test default_cfg.start_inclusive == -1 || cfg.reverse_end_inclusive == -1
    @test default_cfg.end_inclusive == -1 || cfg.reverse_end_inclusive == -1
  end

  # Test manually specified config
  @testset "hard-coded and loaded config" begin
    hardcoded_cfg = example_config()
    @test typeof(hardcoded_cfg) == Configuration
    loaded_cfg = load_config_from_json("test_data/test_config.json")
    @test typeof(loaded_cfg) == Configuration
    @testset "comparing value of $field" for field in fieldnames(typeof(hardcoded_cfg))
      @test getfield(hardcoded_cfg, field) == getfield(loaded_cfg, field)
    end
  end
end
