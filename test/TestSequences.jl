using Test

using BioSequences

function getSequences(filename)::Vector{FASTQ.Record}
    seqs = Vector{FASTQ.Record}()
    iterator = open(FASTQ.Reader, filename)
    for seq in iterator
        push!(seqs, seq)
    end
    return seqs
end

gattaca_seq = "GATTACAGATTACAactggtGATTACAAAAATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGG"
reverse_complement = "CCCTGGCAAATTCATTTCTTCTAATACTGTATCATCTGCTCCTGTATCTAATAGAGCTTCCTTTAGTTGCCCCCCTATTTTTGTAATCaccagtTGTAATCTGTAATC"
test_slices = [((1, -1, -1, 1, false), gattaca_seq, "Full Seq"),
               ((8, -1, 14, -1, false), "GATTACA", "Normal"),
               ((14, -1, 8, -1, false), "ACATTAG", "Backwards"),
               ((-1, 6, -1, 1, false), "CCAGGG", "End bit"),
               ((-1, 4, -1, 8, false), "ACCGT", "End bit backwards"),
               ((1, -1, -1, 1, true), reverse_complement , "Reverse Complement"),
               ((7, -1, 11, -1, true), "CAAAT", "Reverse Complement Slice"),
]

@testset "Fastq Sequences" begin
    seqs = getSequences("test_data/basic.fastq")
    @test length(seqs) == 2
    @test string(FASTQ.sequence(seqs[1])) == uppercase(gattaca_seq)
    @testset "Slicing: $(t[3])" for t in test_slices
        sliced_sequence, sliced_quality = slice_sequence(seqs[1], t[1]...)
        @test string(sliced_sequence) == uppercase(t[2])
        # TODO: Test that the qualitie's get sliced correctly too
    end

end
