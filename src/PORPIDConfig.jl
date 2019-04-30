push!(LOAD_PATH, ".")
module PORPIDConfig
export Configuration, Template, fasta, fastq, read_from_json
import JSON

using States

struct Template
  name::AbstractString
  reference::Vector{AbstractState}
end

Template(name::AbstractString, reference::AbstractString) =
    Template(name, string_to_state_array(reference))

@enum FileType fastq=1 fasta=2

function Base.convert(::Type{FileType}, str_repr::AbstractString)
  return (lowercase(str_repr) == "fasta") ? fasta : fastq
end

mutable struct Configuration
  files::Vector{AbstractString}
  filetype::FileType
  start_inclusive::Integer
  reverse_start_inclusive::Integer
  end_inclusive::Integer
  reverse_end_inclusive::Integer
  max_allowed_errors::Integer
  try_reverse_complement::Bool
  templates::Vector{Template}
end

Configuration() = Configuration(Vector{AbstractString}(undef, 0), fastq, -1, -1, -1, -1, 4, false, Vector{Template}(undef, 0))

function read_from_json(json_file_location)
  config = Configuration()
  params = JSON.parsefile(json_file_location)
  config.files = params["files"]
  if haskey(params, "options")
    config.try_reverse_complement = get(params["options"], "do_reverse_complement", config.try_reverse_complement)
  end
  config.filetype = get(params, "filetype", config.filetype)
  config.start_inclusive = get(params, "start_inclusive", config.start_inclusive)
  config.reverse_start_inclusive = get(params, "reverse_start_inclusive", config.reverse_start_inclusive)
  config.end_inclusive = get(params, "end_inclusive", config.end_inclusive)
  config.reverse_end_inclusive = get(params, "reverse_end_inclusive", config.reverse_end_inclusive)
  config.max_allowed_errors = get(params, "max_allowed_errors", config.max_allowed_errors)
  for json_template in params["templates"]
    template = Template(json_template["name"], json_template["reference"])
    push!(config.templates, template)
  end
  return config
end

end
