push!(LOAD_PATH, ".")
module HMMIDMethods
export process, process_file, sequence_to_observation, best_of_forward_and_reverse, slice_sequence

using NeedlemanWunsch
using HMMIDConfig
using States
using BioSequences
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
function process_file(file_name, config, output_function; print_every=0, print_callback=x->println("Processed $(x) sequences"))
  start_i = config.start_inclusive
  r_start_i = config.reverse_start_inclusive
  end_i = config.end_inclusive
  r_end_i = config.reverse_end_inclusive
  try_reverse_complement = config.try_reverse_complement
  if config.filetype == fastq
    iterator = open(FASTQ.Reader, file_name)
  else
    iterator = open(FASTA.Reader, file_name)
  end
  i = 0
  for sequence in iterator
    if typeof(sequence) <: FASTA.Record
      # Could do this is the slice_sequence method, but if we want to print out the sequences, it's best to have it here
      sequence = FASTQ.Record(FASTQ.identifier(sequence),
                              FASTQ.description(sequence), 
                              FASTQ.sequence(sequence), 
                              fill(Int8(DEFAULT_QUALITY), length(FASTQ.sequence(sequence))))
    end
    forward_seq, forward_quality = slice_sequence(sequence, start_i, r_start_i, end_i, r_end_i, false)
    is_reverse_complement = false
    if try_reverse_complement
      reverse_seq, reverse_quality = slice_sequence(sequence, start_i, r_start_i, end_i, r_end_i, true)
      best_score, best_template, best_tag, best_errors, is_reverse_complement = best_of_forward_and_reverse(forward_seq, forward_quality, reverse_seq, reverse_quality, config.templates)
    else
      best_score, best_template, best_tag, best_errors = choose_best_template(forward_seq, forward_quality, config.templates)
    end
    if is_reverse_complement
      sequence.data[sequence.sequence] = Vector{UInt8}(String(reverse_complement(FASTQ.sequence(sequence))))
      sequence.data[sequence.quality] = reverse(sequence.data[sequence.quality])
    end
    tag = length(best_tag) > 0 ? string(best_tag) : "NO_TAG"
    tag = best_errors <= config.max_allowed_errors ? tag : "REJECTS"
    output_function(file_name, best_template, tag, sequence, best_score)
    i += 1
    if print_every > 0 && i % print_every == 0
      print_callback(i)
    end
  end
end

function slice_sequence(sequence, start_i, r_start_i, end_i, r_end_i, do_reverse_complement)
  if start_i < 0 && r_start_i > 0
    start_i = length(FASTQ.sequence(sequence)) - r_start_i
  end
  if end_i < 0 && r_end_i > 0
    end_i = length(FASTQ.sequence(sequence)) - r_end_i
  end
  start_i = min(length(FASTQ.sequence(sequence)), max(1, start_i))
  end_i = min(length(FASTQ.sequence(sequence)), max(1, end_i))
  # If the start and end are the other way around, we want our sequence to go backwards
  reverse_seq = false
  if start_i > end_i
    reverse_seq = true
    start_i, end_i = end_i, start_i
  end
  if do_reverse_complement
    seq = reverse_complement(FASTQ.sequence(sequence))[start_i:end_i]
    quality = view(sequence.quality, length(sequence.quality)-start_i+1:-1:length(sequence.quality)-end_i+1)
  else
    seq = FASTQ.sequence(sequence)[start_i:end_i]
    quality = view(FASTQ.quality(sequence), start_i:end_i)
  end
  if reverse_seq
    seq = reverse(seq)
    quality = view(quality, length(quality):-1:1)
  end
  return seq, quality
end

#start_i, end_i -> -end_i, -start_i
function tail_indices(start_index, end_index, length)
  return (py_index_to_julia(-end_index, length, true), py_index_to_julia(-start_index, length, true))
end

function best_of_forward_and_reverse(forward_seq, forward_quality, reverse_seq, reverse_quality, templates)
  forward_best_score, forward_best_template, forward_best_tag, forward_best_errors = choose_best_template(forward_seq, forward_quality, templates)
  reverse_best_score, reverse_best_template, reverse_best_tag, reverse_best_errors = choose_best_template(reverse_seq, reverse_quality, templates)
  if forward_best_score > reverse_best_score
    return forward_best_score, forward_best_template, forward_best_tag, forward_best_errors, false
  else
    return reverse_best_score, reverse_best_template, reverse_best_tag, reverse_best_errors, true
  end
end

function choose_best_template(seq, quality, templates)
  best_score = -Inf
  best_template = nothing
  best_tag = nothing
  best_errors = nothing
  for template in templates
    score, tag, errors = NeedlemanWunsch.extract_tag(seq, quality, template.reference)
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
  writer = FASTQ.Writer(fo)
  write(writer, output_sequence)
  close(writer)
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
