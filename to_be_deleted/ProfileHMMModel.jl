push!(LOAD_PATH, ".")
module ProfileHMMModel

using Nucleotides
using Observations
using States

export extract_tag

const INSERTION_START_PROB = 0.01
const INSERTION_EXTEND_PROB = 0.1
const DELETION_START_PROB = 0.01
const DELETION_EXTEND_PROB = 0.1

const L_PROBABILITY_OF_INSERTION = log(PROBABILITY_OF_INSERTION)
const L_PROBABILITY_OF_DELETION = log(PROBABILITY_OF_DELETION)
const L_PROBABILITY_OF_NORMAL_TRANSITION = log(1 - PROBABILITY_OF_DELETION - PROBABILITY_OF_INSERTION)
const L_PROBABILITY_PER_EXTRA_BASE = log(0.05)
const L_NORMAL_TRANSITION = log(1 - (INSERTION_START_PROB + DELETION_START_PROB))
const L_START_INSERTION = log(INSERTION_START_PROB)
const L_EXTEND_INSERTION = log(INSERTION_EXTEND_PROB)
const L_END_INSERTION = log(1 - INSERTION_EXTEND_PROB)
const L_START_DELETION = log(DELETION_START_PROB)
const L_EXTEND_DELETION = log(DELETION_EXTEND_PROB)
const L_END_DELETION = log(1 - DELETION_EXTEND_PROB)

function mismatch(a::DNASymbol, b::DNASymbol)
  if length(intersect(constituents(a), constituents(b))) > 0
    return 0
  else
    return 1
  end
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

immutable TagTrace
  score::Float64
  errors::Int64
  tag::Array{DNASymbol,1}
end

TagTrace(score::Float64) = TagTrace(score, 0, Array{DNASymbol,1}())
TagTrace(score::Float64, errors::Int64) = TagTrace(score, errors, Array{DNASymbol,1}())

function Base.(:+)(tt::TagTrace, score::Float64)
  TagTrace(tt.score + score, tt.errors, tt.tag)
end

function Base.(:+)(tt::TagTrace, score_error::Tuple{Float64, Int64})
  TagTrace(tt.score + score_error[1], tt.errors + score_error[2], tt.tag)
end

function Base.(:+)(tt::TagTrace, symbol::DNASymbol)
  array = copy(tt.tag)
  Base.push!(array, symbol)
  TagTrace(tt.score, tt.errors, array)
end

Base.isless(tt1::TagTrace, tt2::TagTrace) = Base.isless(tt1.score, tt2.score)

  # "nnnnnnnnnnCAGTTTAACTTTTGGGCCATCCATTCC*"
  #         n   n    n    n    n    n             n    n
  #      - I1  I2   I3   I4   I5   I6            I9   I10
  #    /    \ /  \ /  \ /  \ /  \ /  \    n     /  \ /  \
  # Start - S1 - S2 - S3 - S4 - S5 - S6 - S* - S8 - S9 - End
  #     \      X    X    X    X    X    X    X    X    /
  #       - D1 - D2 - D3 - D4 - D5 - D6 - D* - D8 - D9


function extract_tag(observations::Array{Observation,1}, expected::Array{AbstractState, 1})
  # Probability of most likely path to end at state i, on observation j
  # S[i, j] = max(S[i-1, j-1] * normal_transition^, I[i, j-1] * end_insertion, D[i-1, j] * end_deletion^)
  # I[i, j] = max(S[i-1, j-1] * start_insertion, I[i, j-1] * continue_inserting)
  # D[i, j] = max(S[i-1, j-1] * start_of_deletion, D[i-1, j] * continue_deleting)
  # S*[i, j] = max(S[i, j-1], S[i-1, j-1], D[i-1, j])
  # D*[i, j] = max(s[i-1, j-1], D[i-1, j])

  expected = sub(expected, 2:length(expected))
  n = length(expected)

  D = fill(TagTrace(-Inf), n) # Start with log(0) = -Inf probability
  S = fill(TagTrace(-Inf), n)
  I = fill(TagTrace(-Inf), n+1)

  # Values for first observation
  D[1] = TagTrace(L_START_DELETION, 1)
  for i in 2:length(D)
    D[i] = D[i-1] + (L_EXTEND_DELETION, 1)
  end
  S[1] = TagTrace(L_NORMAL_TRANSITION + log(prob(expected[1].value, observations[1].value, observations[1].prob)),
                  mismatch(expected[1].value, observations[1].value))
  for i in 2:length(S)
    if typeof(expected[i]) <: AbstractRepeatingAnyState
      S[i] = D[i-1]
    else
      S[i] = D[i-1] + L_END_DELETION +
        (log(prob(expected[i].value, observations[1].value, observations[1].prob)),
         mismatch(expected[i].value, observations[1].value))
    end
  end
  I[1] = TagTrace(L_START_INSERTION + log(prob(Nucleotides.DNA_N, observations[1].value, observations[1].prob)), 1)
  for i in 1:length(S)
    if typeof(expected[i]) <: AbstractBarcodeState
      S[i] += observations[1].value
      I[i] += observations[1].value
    elseif typeof(expected[i-1]) <: AbstractBarcodeState
      I[i] += observations[1].value
    end
  end

  # Thereafter...

  D[1] = TagTrace(-Inf) # After first observation
  for observation in observations[2:end]
    for i in 2:length(D)
      if typeof(expected[i]) <: AbstractRepeatingAnyState
        D[i] = max(S[i-1], D[i-1]) # A repeating any state can be skipped with no penalty
      else
        D[i] = max(S[i-1] + (L_START_DELETION, 1), D[i-1] + (L_EXTEND_DELETION, 1))
      end
    end

    I[n+1] = max(I[n+1] + (L_EXTEND_INSERTION, 1), S[n] + (L_START_INSERTION, 1))
    for i in length(S):-1:2 # Backwards because we want to use the previous versions of s and i
      if typeof(expected[i]) <: AbstractRepeatingAnyState
        S[i] = max(S[i], S[i-1], D[i-1]) +
          log(prob(Nucleotides.DNA_N, observation.value, observation.prob))
        # No Insertion
      elseif typeof(expected[i-1]) <: AbstractRepeatingAnyState
        S[i] = max(S[i-1], D[i-1]) +
               (log(prob(expected[i].value, observation.value, observation.prob)),
                mismatch(expected[i].value, observation.value))
        # No insertion
      else
        S[i] = max(S[i-1] + L_NORMAL_TRANSITION, I[i] + L_END_INSERTION, D[i-1] + L_END_DELETION) +
          (log(prob(expected[i].value, observation.value, observation.prob)),
           mismatch(expected[i].value, observation.value))
        I[i] = max(I[i] + L_EXTEND_INSERTION, S[i-1] + L_START_INSERTION) +
          (log(prob(Nucleotides.DNA_N, observation.value, observation.prob)),
           1)
      end
    end
    S[1] = I[1] + (L_END_INSERTION + log(prob(expected[1].value, observation.value, observation.prob)),
      mismatch(expected[1].value, observation.value))
    I[1] = I[1] + (L_EXTEND_INSERTION,
      1)
    for i in 1:length(S)
      if typeof(expected[i]) <: AbstractBarcodeState
        S[i] += observation.value
        I[i] += observation.value
      elseif typeof(expected[i-1]) <: AbstractBarcodeState
        I[i] += observation.value
      end
    end
  end
  best_tag_trace = max(I[n+1], S[n], D[n])
  return (best_tag_trace.score, best_tag_trace.tag, best_tag_trace.errors)
end
end
