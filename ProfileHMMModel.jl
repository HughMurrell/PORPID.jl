push!(LOAD_PATH, ".")
module ProfileHMMModel

using Nucleotides
using Observations

export string_to_state_model, viterbi

const INSERTION_START_PROB = 0.01
const INSERTION_EXTEND_PROB = 0.1
const DELETION_START_PROB = 0.01
const DELETION_EXTEND_PROB = 0.1

const L_PROBABILITY_OF_INSERTION = log(PROBABILITY_OF_INSERTION)
const L_PROBABILITY_OF_DELETION = log(PROBABILITY_OF_DELETION)
const L_PROBABILITY_OF_NORMAL_TRANSITION = log(1 - PROBABILITY_OF_DELETION - PROBABILITY_OF_INSERTION)
const L_PROBABILITY_PER_EXTRA_BASE = log(0.05)
const l_normal_transition = log(1 - (insertion_start_prob + deletion_start_prob))
const l_start_insertion = log(insertion_start_prob)
const l_extend_insertion = log(insertion_extend_prob)
const l_end_insertion = log(1 - insertion_extend_prob)
const l_start_deletion = log(deletion_start_prob)
const l_extend_deletion = log(deletion_extend_prob)
const l_end_deletion = log(1 - deletion_extend_prob)

type ObservableState
  index::Int64
  value::DNASymbol
end

type RepeatingAnyState
  index::Int64
end

type ProfileStateModel
  n::Int64
  expected::Array{States, 1}
end

function string_to_state_model(string_sequence::AbstractString)
  # "NNNNNNNNNNCAGTTTAACTTTTGGGCCATCCATTCC*"
  #         n   n    n    n    n    n             n    n
  #      - I1  I2   I3   I4   I5   I6            I9   I10
  #    /    \ /  \ /  \ /  \ /  \ /  \    n     /  \ /  \
  # Start - S1 - S2 - S3 - S4 - S5 - S6 - S* - S8 - S9 - End
  #     \      X    X    X    X    X    X    X    X    /
  #       - D1 - D2 - D3 - D4 - D5 - D6 - D* - D8 - D9

  return ProfileStateModel(length(string_sequence), string_sequence)
end

function getDefault(array::Array{Float64,1}, index::Int64, default::Float64)
  if index < 1
    return default
  else
    return array[index]
  end
end

function get0(array::Array{Float64,1}, index::Int64)
  return getDefault(array, index, 0.0)
end

type TagTrace
  score::Float64
  tag::Array{DNASymbol,1}
end

TagTrace(score::Float64) = TagTrace(score, Array{DNASymbol,1}())

function +(tt::TagTrace, float::Float64)
  TagTrace(tt.score + float, tt.tag)
end

function push!(tt::TagTrace, symbol::DNASymbol)
  array = copy(tt.tag)
  push!(array, symbol)
  tt.tag = array
end

Base.isless(tt1::TagTrace, tt2::TagTrace) = Base.isless(tt1.score, tt2.score)

function extract_tag(observations::Array{Observation,1}, model::ProfileStateModel)
  # Probability of most likely path to end at state i, on observation j
  # S[i, j] = max(S[i-1, j-1] * normal_transition^, I[i, j-1] * end_insertion, D[i-1, j] * end_deletion^)
  # I[i, j] = max(S[i-1, j-1] * start_insertion, I[i, j-1] * continue_inserting)
  # D[i, j] = max(S[i-1, j-1] * start_of_deletion, D[i-1, j] * continue_deleting)
  # S*[i, j] = max(S[i, j-1], S[i-1, j-1], D[i-1, j])
  # D*[i, j] = max(s[i-1, j-1], D[i-1, j])

  D = fill(TagTrace(-Inf), n) # Start with log(0) = -Inf probability
  S = fill(TagTrace(-Inf), n)
  I = fill(TagTrace(-Inf), n+1)

  # Values for first observation
  D[1] = TagTrace(l_start_deletion)
  for i in 2:length(D)
    D[i] = D[i-1] + l_extend_deletion
  end
  S[1] = TagTrace(l_normal_transition + prob(model.states[1].value, observations[1].value, observations[1].prob))
  for i in 2:length(S)
    S[i] = D[i-1] + l_end_deletion + prob(model.states[i].value, observations[1].value, observations[1].prob)
  end
  I[1] = TagTrace(l_start_insertion + prob(DNA_N, observations[1].value, observations[1].prob))

  # Thereafter...

  D[1] = -Inf # After first observation
  for observation in observations[2:end]
    for i in 2:length(D)
      if model.states[i] <: RepeatingAnyState
        D[i] = max(S[i-1], D[i-1]) # A repeating any state can be skipped with no penalty
      else
        D[i] = max(S[i-1] + l_start_deletion, D[i-1] + l_extend_deletion)
      end
    end

    I[n+1] = max(I[n+1] + l_extend_insertion, S[n] + l_start_insertion)
    for i in length(S):-1:2 # Backwards because we want to use the previous versions of s and i
      if model.states[i] <: RepeatingAnyState
        S[i] = max(S[i], S[i-1], D[i-1]) +
          prob(DNA_N, observation.value, observation.prob)
        # No Insertion
      elseif model.states[i-1] <: RepeatingAnyState
        S[i] = max(S[i-1], D[i-1]) +
               prob(model.states[i].value, observation.value, observation.prob)
        # No insertion
      else
        S[i] = max(S[i-1] + l_normal_transition, I[i] + l_end_insertion, D[i-1] + l_end_deletion) +
          prob(model.states[i].value, observation.value, observation.prob)
        I[i] = max(I[i] + l_extend_insertion, S[i-1] + l_start_insertion) +
          prob(DNA_N, observation.value, observation.prob)
      end
      if model.states[i] <: BarcodeState
        push!(S[i], observation.value
      end
    end
    S[1] = I[1] + l_end_insertion
    I[1] = I[1] + l_extend_insertion
  end
  best_tag_trace = max(I[n+1], S[n], D[n])
  return (best_tag_trace.score, best_tag_trace.tag)
end
