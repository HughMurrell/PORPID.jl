push!(LOAD_PATH, ".")

module Observations
using Bio
using Bio.Seq
export Observation, sequence_to_observations, phred_score_to_prob, prob
export PROBABILITY_OF_INSERTION, PROBABILITY_OF_DELETION, L_PROBABILITY_OF_INSERTION, L_PROBABILITY_OF_DELETION
export L_PROBABILITY_PER_EXTRA_BASE, L_PROBABILITY_OF_NORMAL_TRANSITION

const PROBABILITY_OF_INSERTION = 0.01
const PROBABILITY_OF_DELETION = 0.01
const L_PROBABILITY_OF_INSERTION = log(PROBABILITY_OF_INSERTION)
const L_PROBABILITY_OF_DELETION = log(PROBABILITY_OF_DELETION)
const L_PROBABILITY_OF_NORMAL_TRANSITION = log(1 - PROBABILITY_OF_DELETION - PROBABILITY_OF_INSERTION)
const L_PROBABILITY_PER_EXTRA_BASE = log(0.05)

function phred_score_to_prob(score)
  # According to Wikipedia, score = -10*log10(e) where e is probability of base being wrong
  prob_wrong = exp10(score * -0.1)
  return 1 - prob_wrong
end

function fancy_prob(a::DNANucleotide, probA::Float64, b::DNANucleotide, probB::Float64)
  a_b, a_nb, na_b, na_nb = 0, 0, 0, 0
  A = Bio.Seq.compatbits(a)
  B = Bio.Seq.compatbits(b)

  mask = 0x01
  while mask <= 0x08
    a_b += (A & B & mask != 0)
    a_nb += (A & ~B & mask != 0)
    na_b += (~A & B & mask != 0)
    na_nb += (~A & ~B & mask != 0)
    mask <<= 1
  end

  norm_a = probA/(a_b + a_nb)
  norm_na = (1 - probA)/(na_b + na_nb)
  if !isfinite(norm_na)
    norm_na = 0.0
  end
  norm_b = probB/(a_b + na_b)
  norm_nb = (1 - probB)/(a_nb + na_nb)
  if !isfinite(norm_nb)
    norm_nb = 0.0
  end

  return a_b * norm_a * norm_b + a_nb * norm_a * norm_nb +
         na_b * norm_na * norm_b + na_nb * norm_na * norm_nb
end

function prob(expected::DNANucleotide, prob_expected::Float64, observed::DNANucleotide, prob_observed::Float64)
  return fancy_prob(expected, prob_expected, observed, prob_observed)
end

function prob(expected::DNANucleotide, observed::DNANucleotide, prob_observed::Float64)
  return prob(expected, 1.0, observed, prob_observed)
end

end
