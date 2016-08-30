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

      str_tag = join(map(string, best_tag), "")
      if length(str_tag) == 0
        str_tag = "NO_TAG"
      end
      if output_to_file
        str_seq = join(map(string, sequence.seq), "")
        str_quality = join(map(quality_to_char, sequence.quality), "")
        folder_name = splitext(basename(file_name))[1]
        
        if !isdir("$OUTPUT_FOLDER")
          mkdir("$OUTPUT_FOLDER")
        end
        if !isdir("$OUTPUT_FOLDER/$folder_name")
          mkdir("$OUTPUT_FOLDER/$folder_name")
        end
        if !isdir("$OUTPUT_FOLDER/$folder_name/$best_plex_name")
          mkdir("$OUTPUT_FOLDER/$folder_name/$best_plex_name")
        end

        if best_plex_score > best_plex["tag_length"] * log(0.25) + SCORE_THRESHOLD
          open("$OUTPUT_FOLDER/$folder_name/$best_plex_name/$str_tag.fastq", "a") do f
            write(f, "@$(sequence.label)($(round(best_plex_score, 2)))\n$str_seq\n+$(sequence.label)\n$str_quality\n")
          end
        else
          open("$OUTPUT_FOLDER/$folder_name/$best_plex_name/REJECTS.fastq", "a") do f
            write(f, "@$(sequence.label)($(round(best_plex_score, 2)))\n$str_seq\n+$(sequence.label)\n$str_quality\n")
          end
        end
      end

      printif(params, "print_plex", "\t$(best_plex_name)")
      printif(params, "print_tag", "\t$str_tag")
      printif(params, "print_score", "\t$(round(best_plex_score, 2))")
      println()
    end
  end
end

process(ARGS[1])
end
