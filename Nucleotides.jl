import Bio
import Bio.Seq
import Bio.Seq.DNANucleotide

@enum DNANucleotide DNA_A=1 DNA_C=2 DNA_G=3 DNA_T=4
@enum DNANucCombo DNA_R=1 DNA_Y=2 DNA_M=3 DNA_K=4 DNA_S=5 DNA_W=6 DNA_H=7 DNA_B=8 DNA_V=9 DNA_D=10 DNA_N=11

letter_to_nucleotide = Dict([(repr(DNANucleotide(i))[5], DNANucleotide(i)) for i in 1:4])
letter_to_nuc_combo = Dict([(repr(DNANucCombo(i))[5], DNANucCombo(i)) for i in 1:11])

DNASymbol = Union{DNANucleotide, DNANucCombo}

function prob(expected::DNASymbol, prob_expected::Float64, observed::DNASymbol, prob_observed::Float64)
  norm_expected = prob_expected / length(constituents(expected))
  norm_not_expected = (1 - prob_expected) / (4 - length(constituents(expected)))
  norm_observed = prob_observed / length(constituents(observed))
  norm_not_observed = (1 - prob_observed) / (4 - length(constituents(observed)))
  total_prob = 0.0
  for nuc in [DNA_A, DNA_C, DNA_G, DNA_T]
    total_prob += (nuc in constituents(expected) ? norm_expected : norm_not_expected) *
                  (nuc in constituents(observed) ? norm_observed : norm_not_observed)
  end
  return total_prob
end

function prob(expected::DNASymbol, observed::DNASymbol, prob_observed::Float64)
  return prob(expected, 1.0, observed, prob_observed)
end

function prob(expected::DNASymbol, observed::DNANucleotide, prob_observed::Float64)
  if observed in constituents(expected)
    l = length(constituents(expected))
    return prob_observed / l + ((1 - prob_observed) * (l - 1)) / (3 * l)
  else
    return (1 - prob_observed) / 3
  end
end

function constituents(nucleotide::DNANucleotide)
  return [nucleotide]
end

function constituents(combo::DNANucCombo)
  if combo == DNA_R
    return Set([DNA_A, DNA_G])
  elseif combo == DNA_Y
    return Set([DNA_C, DNA_T])
  elseif combo == DNA_M
    return Set([DNA_A, DNA_C])
  elseif combo == DNA_K
    return Set([DNA_G, DNA_T])
  elseif combo == DNA_S
    return Set([DNA_C, DNA_T])
  elseif combo == DNA_W
    return Set([DNA_A, DNA_T])
  elseif combo == DNA_H
    return Set([DNA_A, DNA_C, DNA_T])
  elseif combo == DNA_B
    return Set([DNA_C, DNA_G, DNA_T])
  elseif combo == DNA_V
    return Set([DNA_A, DNA_C, DNA_G])
  elseif combo == DNA_D
    return Set([DNA_A, DNA_G, DNA_T])
  elseif combo == DNA_N
    return Set([DNA_A, DNA_C, DNA_G, DNA_T])
  else
    return Set([])
  end
end

function Base.string(symbol::DNASymbol)
  return repr(symbol)[5]
end

function convert(DNASymbol, symbol::DNASymbol)
  return symbol
end

function convert(DNASymbol, bio_dna_nucleotide::Bio.Seq.DNANucleotide)
  chr = Char(bio_dna_nucleotide)
  return convert(DNASymbol, chr)
end

function convert(DNASymbol, char::Char)
  char = uppercase(char)
  if char in keys(letter_to_nucleotide)
    return letter_to_nucleotide[char]
  elseif char in keys(letter_to_nuc_combo)
    return letter_to_nuc_combo[char]
  else
    throw(DomainError())
  end
end

immutable FastqIterator
  line_iterator::EachLine
end

type Sequence
  label::ASCIIString
  seq::Array{DNASymbol}
  quality::Array{Int8}
end

function char_to_quality(char)
  return Int(char) - 33
end

FastqIterator(file_name::ASCIIString) = FastqIterator(eachline(open(file_name)))
Base.start(fi::FastqIterator) = start(fi.line_iterator)
Base.done(fi::FastqIterator, state) = done(fi.line_iterator, state)

function Base.next(fi::FastqIterator, state)
  label = nothing
  if typeof(state) <: ASCIIString
    label = state
  end
  line = chomp(next(fi.line_iterator, nothing)[1])
  assert(line[1] == '@')
  label = line[2:end]

  line = chomp(next(fi.line_iterator, nothing)[1])
  sequence = Array{DNASymbol,1}()
  for char in line
    push!(sequence, convert(DNASymbol, char))
  end

  line = chomp(next(fi.line_iterator, nothing)[1])
  assert(line[1] == '+')

  line = chomp(next(fi.line_iterator, nothing)[1])
  quality = Array{Int8,1}()
  for char in line
    push!(quality, char_to_quality(char))
  end

  return Sequence(label, sequence, quality), nothing
end
