push!(LOAD_PATH, ".")

module Observations
using Nucleotides
export Observation, sequence_to_observations, prob_to_fastq_score
export PROBABILITY_OF_INSERTION, PROBABILITY_OF_DELETION, L_PROBABILITY_OF_INSERTION, L_PROBABILITY_OF_DELETION
export L_PROBABILITY_PER_EXTRA_BASE, L_PROBABILITY_OF_NORMAL_TRANSITION

const PROBABILITY_OF_INSERTION = 0.01
const PROBABILITY_OF_DELETION = 0.01
const L_PROBABILITY_OF_INSERTION = log(PROBABILITY_OF_INSERTION)
const L_PROBABILITY_OF_DELETION = log(PROBABILITY_OF_DELETION)
const L_PROBABILITY_OF_NORMAL_TRANSITION = log(1 - PROBABILITY_OF_DELETION - PROBABILITY_OF_INSERTION)
const L_PROBABILITY_PER_EXTRA_BASE = log(0.05)

type Observation
  value::DNASymbol
  prob::Float64
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

function reverse_complement(observations::Array{Observation,1})
  rc = Array{Observation, 1}()
  for i in length(observations):-1:1
    push!(rc, Observation(Nucleotides.dna_complement(observations[i].value), observations[i].prob))
  end
  return rc
end

end
