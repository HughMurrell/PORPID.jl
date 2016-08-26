push!(LOAD_PATH, ".")
module HMMIDs
import JSON

using SmithWaterman
using States
using Nucleotides
using Observations
using ProfileHMMModel
using WrongStateModel

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

@enum ExtractionAlgorithm ALG_SMITH_WATERMAN=1 ALG_HMM=2 PROFILE_HMM=3
extraction_functions = Dict(
    ALG_SMITH_WATERMAN => SmithWaterman.extract_tag,
    ALG_HMM => WrongStateModel.extract_tag,
    PROFILE_HMM => ProfileHMMModel.extract_tag)

function process(json_file)
  params = JSON.parsefile(json_file)

  # Options
  extraction_algorithm = ALG_SMITH_WATERMAN
  do_reverse_complement = true
  if haskey(params, "options")
    options = params["options"]
    if haskey(options, "algorithm")
      if lowercase(options["algorithm"]) == "hmm"
        extraction_algorithm = ALG_HMM
      elseif lowercase(options["algorithm"]) == "profile"
        extraction_algorithm = PROFILE_HMM
      end
    end
    if haskey(options, "do_reverse_complement")
      do_reverse_complement = options["do_reverse_complement"]
    end
  end

  # Convert Patterns to StateSequences
  for section in params["sections"]
    for plex in section["multiplexes"]
      reference_state_array = States.string_to_state_array(plex["reference"])
      plex["reference_state_array"] = reference_state_array
    end
  end

  for file_name in params["files"]
    printif(params, "print_filename", "$(file_name)\n")
    for sequence in FastqIterator(file_name)
      printif(params, "print_sequence", "  $(sequence.label)\n")
      for section in params["sections"]
        printif(section, "print_section", "    $(section["name"])\n")
        start_i = py_index_to_julia(get(section, "start_inclusive", 0), length(sequence.seq), true)
        end_i = py_index_to_julia(get(section, "end_inclusive", -1), length(sequence.seq), true)
        printif(section, "print_subsequence", "$(sequence.seq[start_i:end_i])\n")
        observations = sequence_to_observations(sequence.seq[start_i:end_i], sequence.quality[start_i:end_i])

        # Find best matching plex (group in a multiplexed sample)
        best_plex_score = -Inf
        best_plex = "None"
        best_tag = "None"
        for plex in section["multiplexes"]
          extract_tag_fn = extraction_functions[extraction_algorithm]
          score, tag = extract_tag_fn(observations, plex["reference_state_array"])

          if do_reverse_complement
            rc_observations = Observations.reverse_complement(observations)
            rc_score, rc_tag = extract_tag_fn(rc_observations, plex["reference_state_array"])
            if rc_score > score
              score = rc_score
              tag = rc_tag
            end
          end

          printif(section, "print_all_scores", "$(plex["name"]) $(round(score, 2)) $(join(map(string, tag), ""))\n")
          if score > best_plex_score
            best_plex_score = score
            best_plex = plex["name"]
            best_tag = tag
          end
        end
        printif(section, "print_plex", "\t$(best_plex)")
        printif(section, "print_tag", "\t$(join(map(string, best_tag), ""))")
        printif(section, "print_score", "\t$(round(best_plex_score, 2))")
        println()
      end
    end
  end
end

process(ARGS[1])
end
