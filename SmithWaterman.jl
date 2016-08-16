import JSON
push!(LOAD_PATH, ".")

using Nucleotides
using Observations
using WrongStateModel

function SmithWaterman(observations::Array{Observation,1}, states::Array{State,1})
  smb = convert(DNASymbol, 'A')
  tag = [smb, smb, smb]
  return (0, tag, nothing)
end

function printif(dict, key, string)
  if get(dict, key, false)
    print(string)
  end
end

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

function swprocess(json_file)
  params = JSON.parsefile(json_file)

  # Convert Patterns to StateSequences
  for section in params["sections"]
    for plex in section["multiplexes"]
      states = string_to_state_model(plex["reference"])
      plex["states"] = states
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
          score, tag, hmmsequence = SmithWaterman(observations, plex["states"])
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

swprocess(ARGS[1])
