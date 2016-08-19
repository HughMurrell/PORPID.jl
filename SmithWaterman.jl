import JSON
push!(LOAD_PATH, ".")

using Nucleotides
using Observations
using WrongStateModel

@enum AlignOp OP_MATCH=1 OP_DEL=2 OP_INS=3

function SmithWaterman(observations::Array{Observation,1}, states::Array{State,1}, star_insertion_score::Float64=0)
  rows = size(states)[1]
  cols = size(observations)[1] + 1 #first column is for 'before first symbol' position
  scores = zeros(Float64, rows, cols)
  ops = Array{AlignOp}(rows, cols)
  #Left-most column contains all deletions
  for r = 2:rows
    if typeof(states[r]) <: RepeatingAnyState
      scores[r,1] = scores[r-1,1] #Deletion of RepeatingAnyState is free
    else
      scores[r,1] = scores[r-1,1] + L_PROBABILITY_OF_DELETION
    end
    ops[r,1] = OP_DEL
  end
  #Top row contains all insertions
  for c = 2:cols
    scores[1,c] = scores[1,c-1] + L_PROBABILITY_OF_INSERTION
    ops[1,c] = OP_INS
  end
  #Fill in scores for the rest of the matrix
  for r = 2:rows
    for c = 2:cols
      currentstate = states[r]
      if typeof(currentstate) <: RepeatingAnyState
        #RepeatingAnyState represents any number of insertions at star_insertion_score or can be freely deleted
        scores[r,c] = max(scores[r-1,c], scores[r,c-1] + star_insertion_score)
      else
        matchscore = scores[r-1,c-1] + log(prob(currentstate.value, observations[c-1].value, observations[c-1].prob))
        inscore = scores[r,c-1] + L_PROBABILITY_OF_INSERTION
        delscore = scores[r-1,c] + L_PROBABILITY_OF_DELETION
        #The order of operations is significant in the case of multiple paths with the same score
        #In this case, a single path is chosen based on operation order
        if (matchscore >= inscore && matchscore >= delscore)
          scores[r,c] = matchscore
          ops[r,c] = OP_MATCH
        elseif (inscore >= matchscore && inscore >= delscore)
          scores[r,c] = inscore
          ops[r,c] = OP_INS
        else
          scores[r,c] = delscore
          ops[r,c] = OP_DEL
        end
      end
    end
  end
  #Construct tag from score and operation matrices
  tag = Array{DNASymbol, 1}()
  r = rows
  c = cols
  while (r > 1 || c > 1)
    if (ops[r,c] == OP_MATCH)
      if (states[r].value == Nucleotides.DNA_N)
        insert!(tag, 1, observations[c-1].value)
      end
      r = r - 1
      c = c - 1
    elseif (ops[r,c] == OP_INS)
      c = c - 1
    else
      r = r - 1
    end
  end
  return (scores[rows,cols], tag)
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

  # Options
  star_insertion_cost::Float64 = 0
  if (haskey(params, "options"))
    options = params["options"]
    if (haskey(options, "penalise_star") && options["penalise_star"])
      star_insertion_cost = log(0.25)
    end
  end

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
          score, tag = SmithWaterman(observations, plex["states"], star_insertion_cost)
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
