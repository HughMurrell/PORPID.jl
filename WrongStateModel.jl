push!(LOAD_PATH, ".")

module WrongStateModel
using Nucleotides
using Observations
using States

export extract_tag

type IndexedObservableState
  index::Int64
  value::DNASymbol
end

type IndexedStartingState
  index::Int64
end

type IndexedRepeatingAnyState
  index::Int64
end

IndexedState = Union{IndexedObservableState, IndexedStartingState, IndexedRepeatingAnyState}

function states_to_indexed_states(given_states::Array{State,1})
  indexed_states = Array{IndexedState,1}()
  i = 1
  for given in given_states
    if typeof(given) <: StartingState
      new_state = IndexedStartingState(i)
    elseif typeof(given) <: ObservableState
      new_state = IndexedObservableState(i, given.value)
    elseif typeof(given) <: RepeatingAnyState
      new_state = IndexedRepeatingAnyState(i)
    else
      throw(DomainError())
    end
    push!(indexed_states, new_state)
    i += 1
  end
  return indexed_states
end

function emition_prob(observation::Observation, state::IndexedState)
  if typeof(state) <: IndexedStartingState
    return -Inf
  elseif typeof(state) <: IndexedRepeatingAnyState
    return log(0.25)
  else
    return log(prob(state.value, observation.value, observation.prob))
  end
end

function transition_prob(new_state::IndexedState, previous_state::IndexedState)
  # Using linear model
  if typeof(new_state) <: IndexedStartingState
    return -Inf
  elseif (typeof(new_state) <: IndexedRepeatingAnyState) && previous_state == new_state
    return L_PROBABILITY_OF_NORMAL_TRANSITION
  elseif (typeof(new_state) <: IndexedRepeatingAnyState) && new_state.index < previous_state.index
    return -Inf # Cannot go back into a repeating state that you've left
  elseif (typeof(previous_state) <: IndexedRepeatingAnyState) && new_state.index < previous_state.index
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

function extract_tag(observations::Array{Observation,1}, given_states::Array{State,1})
  states = states_to_indexed_states(given_states)

  previous_best_probabilities = Dict{Int64, Float64}(1 => 1.0)
  previous_tags = Dict{Int64, Array{DNASymbol,1}}()
  previous_path = Dict{Int64, Array{Int64,1}}()

  for obs_i in 1:length(observations)
    # Current Probabilities
    current_best_probabilities = Dict{Int64, Float64}()
    # Current Tags
    current_tags = Dict{Int64, Array{DNASymbol,1}}()
    current_path = Dict{Int64, Array{Int64,1}}()
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
      if (typeof(states[current_state_i]) <: IndexedObservableState) && states[current_state_i].value == Nucleotides.DNA_N
          push!(current_tags[current_state_i], observations[obs_i].value)
      end
      current_path[current_state_i] = copy(get(previous_path, best_previous_state_i, Array{Int64,1}()))
      push!(current_path[current_state_i], current_state_i)
    end
    previous_best_probabilities = current_best_probabilities
    previous_tags = current_tags
    previous_path = current_path
  end
  best_tag = ""
  reference_path = []
  best_prob = -Inf
  for final_state_i in keys(previous_best_probabilities)
    if get(previous_best_probabilities, final_state_i, -Inf) > best_prob
      best_tag = previous_tags[final_state_i]
      reference_path = previous_path[final_state_i]
      best_prob = previous_best_probabilities[final_state_i]
    end
  end
  return (best_prob, best_tag)
end
end
