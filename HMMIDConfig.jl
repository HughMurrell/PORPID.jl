push!(LOAD_PATH, ".")
module HMMIDConfig
export Configuration, Template, read_from_json
import JSON

using States
using Nucleotides

type Template
  name::String
  reference::Array{AbstractState,1}
end

Template(name::String, reference::String) =
    Template(name, string_to_state_array(reference))

type Configuration
  files::Array{String,1}
  start_inclusive::Integer
  end_inclusive::Integer
  max_allowed_errors::Integer
  try_reverse::Bool
  templates::Array{Template,1}
end

Configuration() = Configuration(Array{String}(0), 0, -1, 4, false, Array{Template}(0))

function read_from_json(json_file_location)
  config = Configuration()
  params = JSON.parsefile(json_file_location)
  config.files = params["files"]
  if haskey(params, "options")
    config.try_reverse = get(params["options"], "do_reverse_complement", config.try_reverse)
  end
  config.start_inclusive = get(params, "start_inclusive", config.start_inclusive)
  config.end_inclusive = get(params, "end_inclusive", config.end_inclusive)
  config.max_allowed_errors = get(params, "max_allowed_errors", config.max_allowed_errors)
  for json_template in params["templates"]
    template = Template(json_template["name"], json_template["reference"])
    push!(config.templates, template)
  end
  return config
end

end
