using Test
using BioSequences, FASTX
using PORPID
using Random

const EXTRAWORDS = ["GATTACA", "ATTAC", "TAGACAT", "AAACCCTTTGGG"]

function countwords(wordblock::String)
    count = Dict{String, Float32}()
    for word in split(wordblock)
        count[word] = get(count, word, 0.0) + 1.0
    end
    return count
end

function samplewordblock(wordblock::String, keepwords::Vector{String}=[])
    Random.seed!(42)
    allwords = countwords(wordblock)
    samples = Dict{String, Float32}()
    testwords = [Random.randsubseq(collect(keys(allwords)), 0.8)..., EXTRAWORDS..., keepwords...]
    for key in testwords
        samples[key] = get(allwords, key, 0.0)
    end
    return samples, tag_index_mapping(testwords)...
end

function test_neighbors(tag::String, wordblock::String, neighbor_function, keepwords::Vector{String}=Vector{String}())
        neighbor_sample_count, tag_to_index, index_to_tag = samplewordblock(wordblock, keepwords)
        neighbors = neighbor_function(tag, tag_to_index, PORPID.PacBioErrorModel, 0)
        for (tag_index, neighbor_count) in neighbors
            neighbor_sample_count[index_to_tag[tag_index]] -= neighbor_count
        end
        for (tag, value) in neighbor_sample_count
            @test value == 0.0
        end
end

@testset "Testing tag error model" begin

    @testset "String manipulation" begin
        @test PORPID.insert_at("GATTACA", 4, "A") == "GATATACA"
        @test PORPID.replace_at("GATTACA", 2, "C") == "GCTTACA"
        @test PORPID.without("GATTACA", 6) == "GATTAA"
    end

    @testset "TagToIndexMapping" begin
        tags = ["ACAT", "CAAT", "CATA", "CCAT", "GCAT", "CGAT", "CAGT", "CATG"]
        tag_to_index, index_to_tag = tag_index_mapping(tags)
        for tag in tags
            @test typeof(tag_to_index[tag]) == Int32
            @test index_to_tag[tag_to_index[tag]] == tag
        end
    end

    @testset "Insertion Error Neighbors" begin
        tag = "CAT"
        words = """ACAT CAAT CAAT CATA
                   CCAT CCAT CACT CATC
                   TCAT CTAT CATT CATT
                   GCAT CGAT CAGT CATG"""
        keepwords = ["CAAT", "CCAT"]
        test_neighbors(tag, words, PORPID.insertion_neighbors, keepwords)
    end

    @testset "Mutation Error" begin
        tag = "CAT"
        words = """AAT TAT GAT
                   CCT CTT CGT
                   CAA CAC CAG"""
        keepwords = ["CAT"] # We assume the probability of a mutation occurring
        # excludes the probability that the mutation is a self mutation, whatever that might mean
        test_neighbors(tag, words, PORPID.mutation_neighbors, keepwords)
    end

    @testset "Deletion Error" begin
        tag = "CAT"
        words = """AT CT CA"""
        test_neighbors(tag, words, PORPID.deletion_neighbors)
    end
end
