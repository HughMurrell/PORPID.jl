push!(LOAD_PATH, ".")
module Resolving
export tag_index_mapping, prob_observed_tags_given_reals, index_counts

using CustomLDA

const REJECT_TAG = "REJECTS"
const ERROR_RATE = 0.01
const DELETION_RATIO = 0.4
const INSERTION_RATIO = 0.4
const MUTATION_RATIO = 0.2

function process(path)
  path = normpath(path)
  println(STDERR, "Reading tag files...")
  counts = tag_counts(path)
  tag_file_names = tag_to_filename(path)

  println(STDERR, "Generating index mapping...")
  tag_to_index, index_to_tag = tag_index_mapping(Set(keys(counts)))

  println(STDERR, "Generating likelihood distributions...")
  probabilities_array = prob_observed_tags_given_reals(tag_to_index)

  indexed_counts = index_counts(counts, tag_to_index)

  println(STDERR, "Iterating..."  )
  most_likely_real_for_each_obs = CustomLDA.LDA(probabilities_array, indexed_counts)

  tag_count = length(most_likely_real_for_each_obs)
  for observed_index in 1:tag_count
    real_index, prob = most_likely_real_for_each_obs[observed_index]
    observed_tag = index_to_tag[observed_index]
    real_tag = index_to_tag[real_index]
    if (prob > 0.99 && real_index == observed_index)
      println("$(tag_file_names[real_tag])")
    end
  end
end

function index_counts(counts, tag_to_index)
  tags = keys(counts)
  indexed_counts = Array{Int32}(length(tags))
  for tag in tags
    indexed_counts[tag_to_index[tag]] = counts[tag]
  end
  return indexed_counts
end

function tag_index_mapping(tags)
  i = 0
  tag_to_index = Dict{String, Int32}()
  index_to_tag = Dict{Int32, String}()
  for t in tags
    i += 1
    tag_to_index[t] = i
    index_to_tag[i] = t
  end
  return tag_to_index, index_to_tag
end

function tag_to_filename(path)
  file_name_of_tag = Dict{String, String}()
  for file in readdir(path)
    file_name, extension = splitext(file)
    if extension == ".fastq" && file_name != REJECT_TAG
      tag = file_name
      split_file_name = split(file_name, '_')
      if length(split_file_name) == 2
        tag = split_file_name[1]
      end
      file_name_of_tag[tag] = file
    end
  end
  return file_name_of_tag
end

function tag_counts(path)
  counts = Dict{String, Int32}()
  for f in readdir(path)
    fname, extension = splitext(f)
    if extension == ".fastq" && fname != REJECT_TAG
      sp = split(fname, '_')
      if length(sp) == 2
        tag, count = sp
        counts[tag] = parse(Int32, count)
      else
        lines = 0
        open("$path$f") do tagfile
          lines = countlines(tagfile)
        end
        if (lines % 4 != 0)
          println(STDERR, "Warning! $f appears to have incorrect format: expected file to contain 4 lines per sequence, and to end on a new line.")
        end
        counts[fname] = round(lines / 4)
      end
    end
  end
  return counts
end

function prob_observed_tags_given_reals(tag_to_index::Dict{String, Int32})
  prob_observed_tags_given_reals = Array{Array{Tuple{Int32, Float32}}}(length(tag_to_index))
  for observed_tag in keys(tag_to_index)
    observed_index = tag_to_index[observed_tag]
    prob_observed_tags_given_reals[observed_index] = prob_observed_tag_given_reals(observed_tag, tag_to_index)
  end
  return prob_observed_tags_given_reals
end

function prob_observed_tag_given_reals(observed_tag::String, tag_to_index::Dict{String, Int32})
  prob_given_reals_dict = Dict{Int32, Float32}()
  ins_nbrs = insertion_neighbours(observed_tag, tag_to_index)
  for t in ins_nbrs
    prob_given_reals_dict[t] = get(prob_given_reals_dict, t, 0.0) + ERROR_RATE * INSERTION_RATIO * (1/4)
  end
  del_nbrs = deletion_neighbours(observed_tag, tag_to_index)
  for t in del_nbrs
    prob_given_reals_dict[t] = get(prob_given_reals_dict, t, 0.0) + ERROR_RATE * DELETION_RATIO
  end
  mut_nbrs = mutation_neighbours(observed_tag, tag_to_index)
  for t in mut_nbrs
    prob_given_reals_dict[t] = get(prob_given_reals_dict, t, 0.0) + ERROR_RATE * MUTATION_RATIO * (1/3)
  end
  prob_given_reals_dict[tag_to_index[observed_tag]] = (1 - ERROR_RATE) ^ length(observed_tag)
  # Collate dictionary into tuple array
  tuple_array = Array{Tuple{Int32, Float32}}(length(prob_given_reals_dict))
  i = 1
  for (k, v) in prob_given_reals_dict
    tuple_array[i] = (k, v)
    i += 1
  end
  return tuple_array
end

#The three neighbour functions below can have duplicate tags in the returned lists
#This is used to accumulate probabilities when there are multiple error paths between two tags
# e.g. there are three different deletions that turn AAA into AA

#tag -> insertion -> neighbours
function insertion_neighbours(tag::String, tag_to_index::Dict{String, Int32})
  neighbours = Array{Int32}(0)
  word = ""
  for c in "ACTG"
    for i in 1:length(tag) + 1
      word = insert_at(tag, i, c)
      if haskey(tag_to_index, word)
        push!(neighbours, tag_to_index[word])
      end
    end
  end
  return neighbours
end

#tag -> deletion -> neighbours
function deletion_neighbours(tag::String, tag_to_index::Dict{String, Int32})
  neighbours = Array{Int32}(0)
  word = ""
  for i in 1:length(tag)
    word = without(tag, i)
    if haskey(tag_to_index, word)
      push!(neighbours, tag_to_index[word])
    end
  end
  return neighbours
end

#tag -> mutation -> neighbours
function mutation_neighbours(tag::String, tag_to_index::Dict{String, Int32})
  neighbours = Array{Int32}(0)
  word = ""
  for c in "ACTG"
    for i in 1:length(tag)
      word = replace_at(tag, i, c)
      if haskey(tag_to_index, word) && word != tag
        push!(neighbours, tag_to_index[word])
      end
    end
  end
  return neighbours
end

function without(str::String, i)
  return "$(str[1:i-1])$(str[i+1:length(str)])"
end

function insert_at(str::String, i, c)
  return "$(str[1:i-1])$c$(str[i:length(str)])"
end

function replace_at(str::String, i, c)
  return "$(str[1:i-1])$c$(str[i+1:length(str)])"
end

if PROGRAM_FILE == @__FILE__
  process(ARGS[1])
end
end
