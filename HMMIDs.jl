push!(LOAD_PATH, ".")
module HMMIDs
import JSON

using SmithWaterman
using States
using Nucleotides
using Observations
using ProfileHMMModel

const SCORE_THRESHOLD = log(0.01) * 2 - 1e^-10
const OUTPUT_FOLDER = "output"

function py_index_to_julia(py_index, length, bound=false)
  if py_index < 0
    py_index += length
  end
  if bound
    return min(length, max(1, py_index))
  else
    return py_index + 1
  end
end

function printif(dict, key, string)
  if get(dict, key, false)
    print(string)
  end
end

function process(json_file)
  params = JSON.parsefile(json_file)

  # Options
  model = SmithWaterman
  do_reverse_complement = true
  output_to_file = true
  if haskey(params, "options")
    options = params["options"]
    if haskey(options, "algorithm")
      if lowercase(options["algorithm"]) == "profile"
        model = ProfileHMMModel
      end
    end
    if haskey(options, "do_reverse_complement")
      do_reverse_complement = options["do_reverse_complement"]
    end
    if haskey(options, "output_to_file")
      output_to_file = options["output_to_file"]
    end
  end

  # Convert Patterns to StateSequences
  for plex in params["multiplexes"]
    reference_state_array = States.string_to_state_array(plex["reference"])
    plex["reference_state_array"] = reference_state_array
    tag_length = 0
    for state in reference_state_array
      if typeof(state) <: AbstractBarcodeState
        tag_length = tag_length + 1
      end
    end
    plex["tag_length"] = tag_length
  end

  for file_name in params["files"]
    printif(params, "print_filename", "$(file_name)\n")

    #plex name => (tag/cluster => sequences with scores)
    #i.e. folder => (file name => file contents)
    plex_to_cluster_map = Dict{ASCIIString, Dict{ASCIIString, Array{Tuple{Nucleotides.Sequence, Float64}, 1}}}()
    for plex in params["multiplexes"]
      plex_to_cluster_map[plex["name"]] = Dict{ASCIIString, Array{Tuple{Nucleotides.Sequence, Float64}, 1}}()
    end

    for sequence in FastqIterator(file_name)
      printif(params, "print_sequence", "  $(sequence.label)\n")
      start_i = py_index_to_julia(get(params, "start_inclusive", 0), length(sequence.seq), true)
      end_i = py_index_to_julia(get(params, "end_inclusive", -1), length(sequence.seq), true)
      printif(params, "print_subsequence", "$(sequence.seq[start_i:end_i])\n")
      observations = sequence_to_observations(sequence.seq[start_i:end_i], sequence.quality[start_i:end_i])

      # Find best matching plex (group in a multiplexed sample)
      best_plex_score = -Inf
      best_plex = Union{}
      best_plex_name = "None"
      best_tag = "None"
      for plex in params["multiplexes"]
        score, tag = model.extract_tag(observations, plex["reference_state_array"])
        if do_reverse_complement
          rc_observations = Observations.reverse_complement(observations)
          rc_score, rc_tag = model.extract_tag(rc_observations, plex["reference_state_array"])
          if rc_score > score
            score = rc_score
            tag = rc_tag
          end
        end

        printif(params, "print_all_scores", "$(plex["name"]) $(round(score, 2)) $(join(map(string, tag), ""))\n")
        if score > best_plex_score
          best_plex_score = score
          best_plex = plex
          best_plex_name = plex["name"]
          best_tag = tag
        end
      end

      str_tag = length(best_tag) > 0 ? join(map(string, best_tag), "") : "NO_TAG"
      cluster_name = best_plex_score > best_plex["tag_length"] * log(0.25) + SCORE_THRESHOLD ? str_tag : "REJECTS"
      cluster_to_sequences = plex_to_cluster_map[best_plex_name]
      sequence_and_score = (sequence, best_plex_score)
      if !haskey(cluster_to_sequences, cluster_name)
        cluster_to_sequences[cluster_name] = Array{Tuple{Nucleotides.Sequence, Float64}, 1}()
      end
      push!(cluster_to_sequences[cluster_name], sequence_and_score)

      printif(params, "print_plex", "\t$best_plex_name")
      printif(params, "print_tag", "\t$str_tag")
      printif(params, "print_score", "\t$(round(best_plex_score, 2))")
      println()
    end #for each sequence

    if output_to_file
      if !isdir("$OUTPUT_FOLDER")
        mkdir("$OUTPUT_FOLDER")
      end
      folder_name = splitext(basename(file_name))[1]
      if !isdir("$OUTPUT_FOLDER/$folder_name")
        mkdir("$OUTPUT_FOLDER/$folder_name")
      end
      println("Writing output to $(abspath("$OUTPUT_FOLDER/$folder_name")) ...")

      for plex_name in keys(plex_to_cluster_map)
        if !isdir("$OUTPUT_FOLDER/$folder_name/$plex_name")
          mkdir("$OUTPUT_FOLDER/$folder_name/$plex_name")
        end

        cluster_to_sequences = plex_to_cluster_map[plex_name]
        for cluster in keys(cluster_to_sequences)
          open("$OUTPUT_FOLDER/$folder_name/$plex_name/$cluster.fastq", "w") do f
            sequence_and_score_array = cluster_to_sequences[cluster]
            for (sequence, score) in sequence_and_score_array
              str_sequence = join(map(string, sequence.seq), "")
              str_quality = join(map(quality_to_char, sequence.quality), "")
              write(f, "@$(sequence.label)($(round(score, 2)))\n$str_sequence\n+$(sequence.label)\n$str_quality\n")
            end #for each sequence
          end #open file
        end #for each cluster (tag)
      end #for each plex
    end #output to files

  end #for each file to process
end #process function

process(ARGS[1])
end
