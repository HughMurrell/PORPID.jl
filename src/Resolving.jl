push!(LOAD_PATH, ".")
module Resolving
export tag_index_mapping, prob_observed_tags_given_reals, index_counts

using CustomLDA

const REJECT_TAG = "REJECTS"

struct ErrorModel
  error_rate::Float64
  deletion_ratio::Float64
  insertion_ratio::Float64
  mutation_ratio::Float64
end

PacBioErrorModel(error_rate=0.005) = ErrorModel(error_rate, 0.4, 0.4, 0.2)
IlluminaErrorModel(error_rate=0.001) = ErrorModel(error_rate, 0.05, 0.05, 0.9)

function process(path, error_model=PacBioErrorModel())
  path = normpath(path)
  println(STDERR, "Reading tag files...")
  counts = tag_counts(path)
  tag_file_names = tag_to_filename(path)

  println(STDERR, "Generating index mapping...")
  tag_to_index, index_to_tag = tag_index_mapping(Set(keys(counts)))

  println(STDERR, "Generating likelihood distributions...")
  probabilities_array = prob_observed_tags_given_reals(tag_to_index, error_model)

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
  indexed_counts = Vector{Int32}(undef, length(tags))
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

# Traverse a directory, finding files corresponding to tags, returning a dictionary with the tags and their assosciated filepaths
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

function prob_observed_tags_given_reals(tag_to_index::Dict{String, Int32}, error_model::ErrorModel, recurse=0)
  prob_observed_tags_given_reals = Vector{Vector{Tuple{Int32, Float32}}}(undef, length(tag_to_index))
  global memoisation = Dict{Tuple{String, Int64}, Vector{Tuple{Int32,Float32}}}()
  for observed_tag in keys(tag_to_index)
    observed_index = tag_to_index[observed_tag]
    prob_observed_tags_given_reals[observed_index] = prob_observed_tag_given_reals(observed_tag, tag_to_index, error_model, recurse)
  end
  return prob_observed_tags_given_reals
end

function prob_observed_tag_given_reals(observed_tag::String, tag_to_index::Dict{String, Int32}, error_model::ErrorModel, recurse=0)
  global memoisation
  if haskey(memoisation, (observed_tag, recurse))
    return memoisation[(observed_tag, recurse)]
  end
  prob_given_reals_dict = Dict{Int32, Float32}()
  ins_nbrs = insertion_neighbors(observed_tag, tag_to_index, error_model, recurse)
  for (index, prob) in ins_nbrs
    prob_given_reals_dict[index] = get(prob_given_reals_dict, index, 0.0) + error_model.error_rate * error_model.insertion_ratio * (1/4) * prob
  end
  del_nbrs = deletion_neighbors(observed_tag, tag_to_index, error_model, recurse)
  for (index, prob) in del_nbrs
    prob_given_reals_dict[index] = get(prob_given_reals_dict, index, 0.0) + error_model.error_rate * error_model.deletion_ratio * prob
  end
  mut_nbrs = mutation_neighbors(observed_tag, tag_to_index, error_model, recurse)
  for (index, prob) in mut_nbrs
    prob_given_reals_dict[index] = get(prob_given_reals_dict, index, 0.0) + error_model.error_rate * error_model.mutation_ratio * (1/3) * prob
  end
  if haskey(tag_to_index, observed_tag)
    prob_given_reals_dict[tag_to_index[observed_tag]] = (1 - error_model.error_rate) ^ length(observed_tag)
  end
  # Collate dictionary into tuple array
  tuple_array = Vector{Tuple{Int32, Float32}}(undef, length(prob_given_reals_dict))
  i = 1
  for (index, prob) in prob_given_reals_dict
    tuple_array[i] = (index, prob)
    i += 1
  end
  memoisation[(observed_tag, recurse)] = tuple_array
  return tuple_array
end

#The three neighbor functions below can have duplicate tags in the returned lists
#This is used to accumulate probabilities when there are multiple error paths between two tags
# e.g. there are three different deletions that turn AAA into AA

#tag -> insertion -> neighbors
function insertion_neighbors(tag::String, tag_to_index::Dict{String, Int32}, error_model, recurse)
  neighbors = Vector{Tuple{Int32, Float32}}(undef, 0)
  word = ""
  for c in "ACTG"
    for i in 1:length(tag) + 1
      word = insert_at(tag, i, c)
      if recurse > 0
        append!(neighbors, prob_observed_tag_given_reals(word, tag_to_index, error_model, recurse-1))
      elseif haskey(tag_to_index, word)
        push!(neighbors, (tag_to_index[word], 1.0))
      end
    end
  end
  return neighbors
end

#tag -> deletion -> neighbors
function deletion_neighbors(tag::String, tag_to_index::Dict{String, Int32}, error_model, recurse)
  neighbors = Vector{Tuple{Int32, Float32}}(undef, 0)
  word = ""
  for i in 1:length(tag)
    word = without(tag, i)
    if recurse > 0
      append!(neighbors, prob_observed_tag_given_reals(word, tag_to_index, error_model, recurse-1))
    elseif haskey(tag_to_index, word)
      push!(neighbors, (tag_to_index[word], 1.0))
    end
  end
  return neighbors
end

#tag -> mutation -> neighbors
function mutation_neighbors(tag::String, tag_to_index::Dict{String, Int32}, error_model, recurse)
  neighbors = Vector{Tuple{Int32, Float32}}(undef, 0)
  word = ""
  for c in "ACTG"
    for i in 1:length(tag)
      word = replace_at(tag, i, c)
      if recurse > 0
        append!(neighbors, prob_observed_tag_given_reals(word, tag_to_index, error_model, recurse-1))
      elseif haskey(tag_to_index, word) && word != tag
        push!(neighbors, (tag_to_index[word], 1.0))
      end
    end
  end
  return neighbors
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
