push!(LOAD_PATH, ".")
module Testing

using Nucleotides
using States

const SYMBOLS = Array{Nucleotides.DNANucleotide, 1}([Nucleotides.DNA_A, Nucleotides.DNA_C, Nucleotides.DNA_G, Nucleotides.DNA_T])

const TEMPLATE = "AAGGnnnnnnnnnnTTAAGGAATCGACG"
const NUM_SEQUENCES = 300
const NUM_IDS = 15
const SEQUENCE_LENGTH = 100
const OUTPUT_FILE = "testoutput.fastq"

const P_INSERT = 0.012
const P_DELETE = 0.012
const P_MUTATE = 0.007

function copy_with_errors(sequence::Array{Nucleotides.DNANucleotide, 1}, barcode_indices::Set{Int64})
  copy = Array{Nucleotides.DNANucleotide, 1}()
  barcode_errors = 0
  primer_errors = 0
  i = 1
  while i <= length(sequence)
    roll = rand(Float64)
    #Insert
    if roll <= P_INSERT
      push!(copy, rand(SYMBOLS))
      if i in barcode_indices || i - 1 in barcode_indices
        barcode_errors += 1
      elseif i <= length(TEMPLATE)
        primer_errors += 1
      end
    #Delete
    elseif roll <= P_DELETE + P_INSERT
      if i in barcode_indices
        barcode_errors += 1
      elseif i <= length(TEMPLATE)
        primer_errors += 1
      end
      i += 1
    #Mutate
    elseif roll <= P_MUTATE + P_DELETE + P_INSERT
      sym = rand(SYMBOLS)
      while sym == sequence[i]
        sym = rand(SYMBOLS)
      end
      push!(copy, sym)
      if i in barcode_indices
        barcode_errors += 1
      elseif i <= length(TEMPLATE)
        primer_errors += 1
      end
      i += 1
    #No Error
    else
      push!(copy, sequence[i])
      i += 1
    end
  end
  #Chance to insert at the end of the sequence
  closed = false
  while !closed
    roll = rand(Float64)
    if roll <= P_INSERT
      push!(copy, rand(SYMBOLS))
      if length(sequence) in barcode_indices
        barcode_errors += 1
      end
    else
      closed = true
    end
  end
  return copy, barcode_errors, primer_errors
end

function generate_test_sequences()
  template_states = States.string_to_state_array(TEMPLATE)
  primers = Array{Array{Nucleotides.DNANucleotide, 1}, 1}()
  barcode_indices = Set{Int64}()
  for i in 1:length(template_states)
    if typeof(template_states[i]) <: AbstractBarcodeState
      push!(barcode_indices, i)
    end
  end
  for p in 1:NUM_IDS
    primer = Array{Nucleotides.DNANucleotide, 1}()
    for state in template_states
      if typeof(state) <: AbstractBarcodeState
        push!(primer, rand(SYMBOLS))
      elseif typeof(state) <: AbstractObservableState
        push!(primer, state.value)
      end
    end
    push!(primers, primer)
  end
  tail = Array{Nucleotides.DNANucleotide, 1}()
  for c in 1:(SEQUENCE_LENGTH - length(TEMPLATE))
    push!(tail, rand(SYMBOLS))
  end

  open(OUTPUT_FILE, "w") do f
    for i in 1:NUM_SEQUENCES
      id = rand(1:length(primers))
      seq, barcode_errors, primer_errors = copy_with_errors(cat(1, primers[id], tail), barcode_indices)
      quality = ""
      for q in 1:length(seq)
        quality = "$(quality)~"
      end
      write(f, "@$id ($barcode_errors, $primer_errors)\n$(join(map(string, seq), ""))\n+\n$quality\n")
    end
  end
end

generate_test_sequences()

end
