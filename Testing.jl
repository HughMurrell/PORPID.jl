push!(LOAD_PATH, ".")
module Testing

using BioSequences
using States

const SYMBOLS = Array{DNA, 1}([DNA_A, DNA_C, DNA_G, DNA_T])

const TEMPLATES = ["AAGGnnnnnnnnnnTTAAGGAATCGACG", "ACACAGTGTACAGTCTGACGTTGCnnnnnnnnCCACTTGCCACCCATBTTATAGCA"]
const NUM_SEQUENCES = 300
const NUM_IDS = 15
const SEQUENCE_LENGTH = 100
const OUTPUT_FILE = "testoutput.fastq"

const P_INSERT = 0.012
const P_DELETE = 0.012
const P_MUTATE = 0.007

function copy_with_errors(sequence::Array{DNA, 1}, template::Array{States.AbstractState}, barcode_indices::Set{Int64})
  copy = Array{DNA, 1}()
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
      elseif i <= length(template)
        primer_errors += 1
      end
    #Delete
    elseif roll <= P_DELETE + P_INSERT
      if i in barcode_indices
        barcode_errors += 1
      elseif i <= length(template)
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
      elseif i <= length(template)
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
  plexes = map(States.string_to_state_array, TEMPLATES)

  #each template has a set of barcode indices
  barcode_indices = Array{Set{Int64}, 1}()
  for template in plexes
    barcode_index_set = Set{Int64}()
    for i in 1:length(template)
      if typeof(template[i]) <: AbstractBarcodeState
        push!(barcode_index_set, i)
      end
    end
    push!(barcode_indices, barcode_index_set)
  end

  #each template has a separate array of primers
  primers = Array{Array{Array{DNA, 1}, 1}, 1}()
  for i in 1:length(plexes)
    push!(primers, Array{Array{DNA, 1}, 1}())
  end
  for p in 1:NUM_IDS
    template_num = rand(1:length(plexes))
    template = plexes[template_num]
    primer = Array{DNA, 1}()
    for state in template
      if typeof(state) <: AbstractBarcodeState
        push!(primer, rand(SYMBOLS))
      elseif typeof(state) <: AbstractObservableState
        push!(primer, state.value)
      end
    end
    push!(primers[template_num], primer)
  end

  #each template has a separate 'tail', i.e. the sequence following the primer
  tails = Array{Array{DNA, 1}, 1}()
  for template in plexes
    tail = Array{DNA, 1}()
    for c in 1:(SEQUENCE_LENGTH - length(template))
      push!(tail, rand(SYMBOLS))
    end
    push!(tails, tail)
  end

  open(OUTPUT_FILE, "w") do f
    for i in 1:NUM_SEQUENCES
      template_num = rand(1:length(plexes))
      primer_num = rand(1:length(primers[template_num]))
      canonical_sequence = cat(1, primers[template_num][primer_num], tails[template_num])
      output_sequence, barcode_errors, primer_errors = copy_with_errors(canonical_sequence, plexes[template_num], barcode_indices[template_num])
      quality = ""
      for q in 1:length(output_sequence)
        quality = "$(quality)~"
      end
      write(f, "@$template_num.$primer_num ($barcode_errors, $primer_errors)\n$(join(map(string, output_sequence), ""))\n+\n$quality\n")
    end
  end
end

generate_test_sequences()

end
