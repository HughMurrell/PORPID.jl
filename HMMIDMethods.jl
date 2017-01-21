push!(LOAD_PATH, ".")
module HMMIDMethods
export process, process_file, sequence_to_observation, best_of_forward_and_reverse

using NeedlemanWunsch
using HMMIDConfig
using States
using Nucleotides
using Observations
using LRUExample

const DEFAULT_MAX_ERRORS = 2
const DEFAULT_MAX_FILE_DESCRIPTORS = 1024
const OUTPUT_FOLDER = "output"

#Go through each template
#get a score
#Method: remember which is best
#Check the number of errors. If too many, tag is "REJECTS"
#Check the tag length, if zero, tag is "NO_TAG"
#If reversed, reverse the sequence
#Push sequence onto tag "cluster"
#Output all clusters to files.

#It should be an iterator for every Sequence of every file listed, converted into a observation array.
#Then run the "Best Score" method with the sequence and the list of potential template's.
#Then save to file

# Json parameters object/type?
# Patterns

# files
function process(json_file_location)
  # Get configuration from json
  config = HMMIDConfig.read_from_json(json_file_location)
  for file_name in config.files
    process_file(file_name, config)
  end
end

# For every file
function process_file(file_name, config)
  start_i = config.start_inclusive
  end_i = config.end_inclusive
  try_reverse = config.try_reverse
  for sequence in FastqIterator(file_name)
    best_score, best_template, best_tag, best_errors, is_reversed = best_of_forward_and_reverse(
        sequence_to_observations(sequence, start_i, end_i, try_reverse), config.templates)
    if is_reversed
      sequence = reverse_complement(sequence)
    end
    tag = length(best_tag) > 0 ? join(map(string, best_tag), "") : "NO_TAG"
    tag = best_errors <= config.max_allowed_errors ? tag : "REJECTS"
    write_to_file(file_name, tag, sequence, best_score)
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

function write_to_file(source_file_name, tag, output_sequence, score)
  println("$(output_sequence.label) $(tag) $(score)")
  # Get output file name
  # Check LRU for existing reference to file
  # If yes, use reference, if not, create new one and store in LRU
  # Append this sequence to the file
end

end
