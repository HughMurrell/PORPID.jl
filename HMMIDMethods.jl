push!(LOAD_PATH, ".")
module HMMIDMethods
export process, process_file, sequence_to_observation, best_of_forward_and_reverse

using NeedlemanWunsch
using HMMIDConfig
using States
using Nucleotides
using Observations

const DEFAULT_MAX_ERRORS = 2
const DEFAULT_MAX_FILE_DESCRIPTORS = 1024
const OUTPUT_FOLDER = "output"
const DEFAULT_QUALITY = 30

function process(json_file_location, output_function)
  # Get configuration from json
  config = HMMIDConfig.read_from_json(json_file_location)
  for file_name in config.files
    process_file(file_name, config, output_function)
  end
end

# For every file
function process_file(file_name, config, output_function)
  start_i = config.start_inclusive
  end_i = config.end_inclusive
  try_reverse = config.try_reverse
  if config.filetype == fastq
    iterator = FastqIterator(file_name)
  else
    iterator = FastaIterator(file_name)
  end
  for sequence in iterator
    if typeof(sequence) <: FastaSequence
      sequence = FastqSequence(sequence.label, sequence.seq, fill(DEFAULT_QUALITY, length(sequence.seq)))
    end
    best_score, best_template, best_tag, best_errors, is_reversed = best_of_forward_and_reverse(
        sequence_to_observations(sequence, start_i, end_i, try_reverse), config.templates)
    if is_reversed
      sequence = reverse_complement(sequence)
    end
    tag = length(best_tag) > 0 ? join(map(string, best_tag), "") : "NO_TAG"
    tag = best_errors <= config.max_allowed_errors ? tag : "REJECTS"
    output_function(file_name, best_template, tag, sequence, best_score)
  end
end

function sequence_to_observations(sequence, start_i, end_i, do_reverse)
  start_i = py_index_to_julia(start_i, length(sequence.seq), true)
  end_i = py_index_to_julia(end_i, length(sequence.seq), true)
  observations = Observations.sequence_to_observations(sequence.seq[start_i:end_i], sequence.quality[start_i:end_i])
  if !do_reverse
    return observations
  end
  tail_start_i, tail_end_i = tail_indices(start_i, end_i, length(sequence.seq))
  tail_observations = Observations.sequence_to_observations(sequence.seq[tail_start_i:tail_end_i], sequence.quality[tail_start_i:tail_end_i])
  rc_observations = Observations.reverse_complement(tail_observations)
  return observations, rc_observations
end

function py_index_to_julia(py_index, length, bound=false)
  if py_index < 0
    py_index += length
  end
  if bound
    return min(length, max(1, py_index + 1))
  else
    return py_index + 1
  end
end

#start_i, end_i -> -end_i, -start_i
function tail_indices(start_index, end_index, length)
  return (py_index_to_julia(-end_index, length, true), py_index_to_julia(-start_index, length, true))
end

function best_of_forward_and_reverse(i, templates)
  if typeof(i) <: Tuple
    forward_observation, reverse_observation = i
    forward_best_score, forward_best_template, forward_best_tag, forward_best_errors = best_template(forward_observation, templates)
    reverse_best_score, reverse_best_template, reverse_best_tag, reverse_best_errors = best_template(reverse_observation, templates)
    if forward_best_score > reverse_best_score
      return forward_best_score, forward_best_template, forward_best_tag, forward_best_errors, false
    else
      return reverse_best_score, reverse_best_template, reverse_best_tag, reverse_best_errors, true
    end
  elseif typeof(i) <: Array{Observation}
    return best_template(i, templates)..., false
  end
end

function best_template(observations, templates)
  best_score = -Inf
  best_template = nothing
  best_tag = nothing
  best_errors = nothing
  for template in templates
    score, tag, errors = NeedlemanWunsch.extract_tag(observations, template.reference)
    if score >= best_score
      best_score, best_template, best_tag, best_errors = score, template, tag, errors
    end
  end
  return best_score, best_template, best_tag, best_errors
end

function write_to_file(source_file_name, template, tag, output_sequence, score)
  source_file_name = basename(source_file_name)
  output_file_name = "output/$(source_file_name)/$(template.name)/$(tag).fastq"
  mkpath("output/$(source_file_name)/$(template.name)")
  fo = open(output_file_name, "a")
  str_sequence = join(map(string, output_sequence.seq), "")
  str_quality = join(map(quality_to_char, output_sequence.quality), "")
  write(fo, "@$(output_sequence.label)($(round(score, 2)))\n$str_sequence\n+\n$str_quality\n")
  close(fo)
end

function write_to_dictionary(dictionary, source_file_name, template, tag, output_sequence, score)
  directory = "$(source_file_name)/$(template.name)"
  if !haskey(dictionary, directory)
    dictionary[directory] = Dict()
  end
  directory_dict = dictionary[directory]
  if !haskey(directory_dict, tag)
    directory_dict[tag] = []
  end
  push!(directory_dict[tag], (score, output_sequence))
end

function write_to_file_count_to_dict(dictionary, source_file_name, template, tag, output_sequence, score)
  write_to_file(source_file_name, template, tag, output_sequence, score)
  directory = "$(source_file_name)/$(template.name)"
  if !haskey(dictionary, directory)
    dictionary[directory] = Dict()
  end
  directory_dict = dictionary[directory]
  directory_dict[tag] = get(directory_dict, tag, 0) + 1
end

if basename(PROGRAM_FILE) == basename(@__FILE__)
  println("Processing $(ARGS[1])")
  dir_dict = Dict()
  my_output_func(source_file_name, template, tag, output_sequence, score) = write_to_file_count_to_dict(dir_dict, source_file_name, template, tag, output_sequence, score)
  HMMIDMethods.process(ARGS[1], my_output_func)
  for dir in keys(dir_dict)
    println(dir)
    sizes = []
    for key in keys(dir_dict[dir])
      push!(sizes, (dir_dict[dir][key], key))
    end
    sort!(sizes)
    for (size, key) in sizes
      println(key, "\t", size)
    end
  end
end
end
