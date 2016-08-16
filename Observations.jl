push!(LOAD_PATH, ".")

module Observations
using Nucleotides
export Observation, sequence_to_observations

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
end
