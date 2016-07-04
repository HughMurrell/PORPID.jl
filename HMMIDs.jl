module HMMIDs
import JSON
include("Nucleotides.jl")

#const MAX_SHIFT = 30
#const MAX_INDEL = 15
const PROBABILITY_OF_INSERTION = 0.01
const PROBABILITY_OF_DELETION = 0.01
const L_PROBABILITY_OF_INSERTION = log(PROBABILITY_OF_INSERTION)
const L_PROBABILITY_OF_DELETION = log(PROBABILITY_OF_DELETION)
const L_PROBABILITY_OF_NORMAL_TRANSITION = log(1 - PROBABILITY_OF_DELETION - PROBABILITY_OF_INSERTION)
const L_PROBABILITY_PER_EXTRA_BASE = log(0.05)

type ObservableState
  index::Int64
  value::DNASymbol
end

type StartingState
  index::Int64
end

type RepeatingAnyState
  index::Int64
end

State = Union{ObservableState, StartingState, RepeatingAnyState}

type Observation
  value::DNASymbol
  prob::Float64
end

function emition_prob(observation::Observation, state::State)
  if typeof(state) <: StartingState
    return -Inf
  elseif typeof(state) <: RepeatingAnyState
    return log(0.25)
  else
    return log(prob(state.value, observation.value, observation.prob))
  end
end

function transition_prob(new_state::State, previous_state::State)
  # Using linear model
  if typeof(new_state) <: StartingState
    return -Inf
  elseif (typeof(new_state) <: RepeatingAnyState) && previous_state == new_state
    return L_PROBABILITY_OF_NORMAL_TRANSITION
  elseif (typeof(new_state) <: RepeatingAnyState) && new_state.index < previous_state.index
    return -Inf # Cannot go back into a repeating state that you've left
  elseif (typeof(previous_state) <: RepeatingAnyState) && new_state.index < previous_state.index
    return -Inf # Cannot go backwards once you've entered a repeating state
  end
  diff = new_state.index - previous_state.index
  if diff == 1
    return L_PROBABILITY_OF_NORMAL_TRANSITION
  elseif diff < 1
    return L_PROBABILITY_OF_INSERTION + -diff * L_PROBABILITY_PER_EXTRA_BASE
  elseif diff > 1
    return L_PROBABILITY_OF_DELETION + (diff - 2) * L_PROBABILITY_PER_EXTRA_BASE
  end
end

function sliding_window(max_range::Symbol, max_length::Int64, index::Int64)
  if isdefined(max_range)
    offset = eval(max_range)::Int64
    return max(1, index - offset):min(max_length, index + offset)
  else
    return 1:max_length
  end
end

function viterbi(observations::Array{Observation,1}, states::Array{State,1})
  previous_best_probabilities = Dict{Int64, Float64}(1 => 1.0)
  previous_tags = Dict{Int64, Array{DNASymbol,1}}()
  previous_sequence = Dict{Int64, Array{Int64,1}}()

  for obs_i in 1:length(observations)
    # Current Probabilities
    current_best_probabilities = Dict{Int64, Float64}()
    # Current Tags
    current_tags = Dict{Int64, Array{DNASymbol,1}}()
    current_sequence = Dict{Int64, Array{Int64,1}}()
    for current_state_i in sliding_window(:MAX_SHIFT, length(states), obs_i)
      best_previous_state_i = -1
      best_path_prob = -Inf
      for previous_state_i in sliding_window(:MAX_INDEL, length(states), current_state_i)
        previous_state_prob = get(previous_best_probabilities, previous_state_i, -Inf)
        path_prob = previous_state_prob + transition_prob(states[current_state_i], states[previous_state_i])
        if path_prob > best_path_prob
          best_previous_state_i = previous_state_i
          best_path_prob = path_prob
        end
      end
      prob = best_path_prob + emition_prob(observations[obs_i], states[current_state_i])
      current_best_probabilities[current_state_i] = prob
      current_tags[current_state_i] = copy(get(previous_tags, best_previous_state_i, Array{DNASymbol,1}()))
      if (typeof(states[current_state_i]) <: ObservableState) && states[current_state_i].value == DNA_N
          push!(current_tags[current_state_i], observations[obs_i].value)
      end
      current_sequence[current_state_i] = copy(get(previous_sequence, best_previous_state_i, Array{Int64,1}()))
      push!(current_sequence[current_state_i], current_state_i)
    end
    previous_best_probabilities = current_best_probabilities
    previous_tags = current_tags
    previous_sequence = current_sequence
  end
  best_tag = ""
  best_sequence = []
  best_prob = -Inf
  for final_state_i in keys(previous_best_probabilities)
    if get(previous_best_probabilities, final_state_i, -Inf) > best_prob
      best_tag = previous_tags[final_state_i]
      best_sequence = previous_sequence[final_state_i]
      best_prob = previous_best_probabilities[final_state_i]
    end
  end
  return (best_prob, best_tag, best_sequence)
end

function string_to_states(string_sequence)
  states = Array{State,1}()
  push!(states, StartingState(1))
  i = 2
  for char in string_sequence
    if char == '*'
      state = RepeatingAnyState(i)
    else
      state = ObservableState(i, convert(DNASymbol, char))
    end
    push!(states, state)
    i += 1
  end
  return states
end

function fastq_score_to_prob(score)
  # According to Wikipedia, score = -10*log10(e) where e is probability of base being wrong
  prob_wrong = exp10(score * -0.1)
  return 1 - prob_wrong
end

function sequence_to_observations(sequence, quality)
  observations = Array{Observation,1}()
  for (base, score) in zip(sequence, quality)
    push!(observations, Observation(convert(DNASymbol, base), fastq_score_to_prob(score)))
  end
  return observations
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

function printif(dict, key, string)
  if get(dict, key, false)
    print(string)
  end
end

function process(json_file)
  params = JSON.parsefile(json_file)

  # Convert Patterns to StateSequences
  for section in params["sections"]
    for plex in section["multiplexes"]
      states = string_to_states(plex["reference"])
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
          score, tag, hmmsequence = viterbi(observations, plex["states"])
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
