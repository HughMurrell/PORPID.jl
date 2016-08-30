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

function process(json_file)
  params = JSON.parsefile(json_file)

  # Options
  model = SmithWaterman
  do_reverse_complement = true
  if haskey(params, "options")
    options = params["options"]
    if haskey(options, "algorithm")
      if lowercase(options["algorithm"]) == "hmm"
        model = WrongStateModel
      elseif lowercase(options["algorithm"]) == "profile"
        model = ProfileHMMModel
      end
    end
    if haskey(options, "do_reverse_complement")
      do_reverse_complement = options["do_reverse_complement"]
    end
  end

  # Convert Patterns to StateSequences
  for plex in params["multiplexes"]
    reference_state_array = States.string_to_state_array(plex["reference"])
    plex["reference_state_array"] = reference_state_array
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
      best_plex = "None"
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
          best_plex = plex["name"]
          best_tag = tag
        end
      end
      printif(params, "print_plex", "\t$(best_plex)")
      printif(params, "print_tag", "\t$(join(map(string, best_tag), ""))")
      printif(params, "print_score", "\t$(round(best_plex_score, 2))")
      println()
    end
  end
end

process(ARGS[1])
end
