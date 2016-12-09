push!(LOAD_PATH, ".")
module SimpleConsensus
using Nucleotides
using Distributions

function read_fasta_records(filename)
    seqs = Array{Array{DNASymbol}}(0)
    for sequence in FastaIterator(filename)
      push!(seqs, sequence.seq)
    end
    return seqs
end

function consensus(fastaPath)
    seqs = read_fasta_records(fastaPath)
    first_length = length(seqs[1])
    for seq in seqs
      if length(seq) != first_length
        error("Aligned sequences must all be the same length.")
      end
    end
    consensus_seq = [mode([seq[i] for seq in seqs]) for i in 1:first_length]
    filter!(e->e!=Nucleotides.DNA_GAP, consensus_seq)
    return consensus_seq
end

if length(ARGS) > 0
  print(join(map(string, consensus(ARGS[1]))))
end
end
