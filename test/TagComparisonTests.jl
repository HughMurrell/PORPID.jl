include("../src/include_all.jl")
using Test
using BioSequences
using Resolving
using Random

const extra_words = ["GATTACA", "ATTAC", "TAGACAT", "AAACCCTTTGGG"]

function count_words(wordblock::String)
    count = Dict{String, Float32}()
    for word in split(wordblock)
        count[word] = get(count, word, 0.0) + 1.0
    end
    return count
end

function sample_word_block(wordblock::String, keep_words::Vector{String}=[])
    Random.seed!(42)
    all_words = count_words(wordblock)
    samples = Dict{String, Float32}()
    test_words = [Random.randsubseq(collect(keys(all_words)), 0.8)..., extra_words..., keep_words...]
    for key in test_words
        samples[key] = get(all_words, key, 0.0)
    end
    return samples, tag_index_mapping(test_words)...
end

function testNeighbors(tag::String, wordblock::String, neighborFunction, keep_words::Vector{String}=Vector{String}())
        neighbor_sample_count, tag_to_index, index_to_tag = sample_word_block(wordblock, keep_words)
        neighbors = neighborFunction(tag, tag_to_index, Resolving.PacBioErrorModel, 0)
        for (tag_index, neighbor_count) in neighbors
            neighbor_sample_count[index_to_tag[tag_index]] -= neighbor_count
        end
        for (tag, value) in neighbor_sample_count
            @test value == 0.0
        end
end

@testset "Testing tag error model" begin

    @testset "String manipulation" begin
        @test Resolving.insert_at("GATTACA", 4, "A") == "GATATACA"
        @test Resolving.replace_at("GATTACA", 2, "C") == "GCTTACA"
        @test Resolving.without("GATTACA", 6) == "GATTAA"
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
        keep_words = ["CAAT", "CCAT"]
        testNeighbors(tag, words, Resolving.insertion_neighbors, keep_words)
    end

    @testset "Mutation Error" begin
        tag = "CAT"
        words = """AAT TAT GAT
                   CCT CTT CGT
                   CAA CAC CAG"""
        keep_words = ["CAT"] # We assume the probability of a mutation occurring
        # excludes the probability that the mutation is a self mutation, whatever that might mean
        testNeighbors(tag, words, Resolving.mutation_neighbors, keep_words)
    end

    @testset "Deletion Error" begin
        tag = "CAT"
        words = """AT CT CA"""
        testNeighbors(tag, words, Resolving.deletion_neighbors)
    end
end
