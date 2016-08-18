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

type ProfileStateModel
  n::Int64
  template_string::AbstractString
end

function string_to_state_model(string_sequence::AbstractString)
  # "NNNNNNNNNNCAGTTTAACTTTTGGGCCATCCATTCC*"
  #         n   n    n    n    n    n   n
  #      - I1  I2   I3   I4   I5   I6  I7
  #    /    \ /  \ /  \ /  \ /  \ /  \/  \
  # Start - S1 - S2 - S3 - S4 - S5 - S6 - End
  #     \      X    X    X    X    X    /
  #       - D1 - D2 - D3 - D4 - D5 - D6

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


function viterbi(observations::Array{Observation,1}, model::ProfileStateModel)
  # Probability of most likely path to end at state i, on observation j
  # S[i, j] = max(S[i-1, j-1] * normal_transition, I[i, j-1] * end_insertion, D[i-1, j] * end_deletion)
  # I[i, j] = max(S[i-1, j-1] * start_insertion, I[i, j-1] * continue_inserting)
  # D[i, j] = max(S[i-1, j-1] * start_of_deletion, D[i-1, j] * continue_deleting)

  D = fill(-Inf, n) # Zero log probability: log(0) = -Inf
  S = fill(-Inf, n)
  I = fill(-Inf, n+1)

  # Values for first observation
  D[1] = l_start_deletion
  for i in 2:length(D)
    D[i] = D[i-1] + l_extend_deletion
  end
  S[1] = l_normal_transition # + observation_prob
  for i in 2:length(S)
    S[i] = D[i-1] + l_end_deletion # + observation_prob
  end
  I[1] = l_start_insertion # + observation_prob

  # Thereafter...

  D[1] = -Inf # After first observation
  for observation in observations[2:end]
    for i in 2:length(D)
      D[i] = max(S[i-1] + l_start_deletion, D[i-1] + l_extend_deletion)
    end

    I[n+1] = max(I[n+1] + l_extend_insertion, S[n] + l_start_insertion)
    for i in length(S):-1:2 # Backwards because we want to use the previous versions of s and i
      S[i] = max(S[i-1] + l_normal_transition, I[i] + l_end_insertion, D[i-1] + l_end_deletion)
        # + observation_prob
      I[i] = max(I[i] + l_extend_insertion, S[i-1] + l_start_insertion)
        # + observation_prob
    end
    S[1] = I[1] + l_end_insertion
    I[1] = I[1] + l_extend_insertion
  end
  return (best_prob, best_tag, reference_path)
end
