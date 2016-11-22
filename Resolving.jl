push!(LOAD_PATH, ".")
module Resolving

using CustomLDA

const REJECT_TAG = "REJECTS"
const ERROR_RATE = 0.005
const DELETION_RATIO = 0.4
const INSERTION_RATIO = 0.4
const MUTATION_RATIO = 0.2

function process(path)
  path = normpath(path)
  println(STDERR, "Reading tag files...")
  @time counts = tag_counts(path)

  println(STDERR, "Generating likelihood distributions...")
  @time probabilities_dictionary = prob_observed_tags_given_reals(counts)

  println(STDERR, "Generating index mapping...")
  @time tag_to_index, index_to_tag = tag_index_mapping(Set(keys(counts)))

  indexed_counts = index_counts(counts, tag_to_index)

  println(STDERR, "Converting tag mapping to row array...")
  @time probabilities_array = probabilities_to_row_array(probabilities_dictionary, tag_to_index)

  println(STDERR, "Iterating..."  )
  @time most_likely_real_for_each_obs = CustomLDA.LDA(probabilities_array, indexed_counts)

  println("tag,count,best_tag,best_score")
  tag_count = length(most_likely_real_for_each_obs)
  for observed_index in 1:tag_count
    real_index, prob = most_likely_real_for_each_obs[observed_index]
    observed_tag = index_to_tag[observed_index]
    real_tag = index_to_tag[real_index]
    #println("$(index_to_tag[r])\t$(round(posterior[r,maxindex], 4))")
    #if (maxval < 0.99)
    println("$(observed_tag),$(counts[observed_tag]),$(real_tag),$(round(prob, 4))")
    #end
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

function probabilities_to_row_array(probabilities, tag_to_index)
  all_tags = keys(probabilities)
  cl = Array{Array{Tuple{Int32, Float32}}}(length(all_tags))
  for observation in all_tags
    observation_index = tag_to_index[observation]
    cl[observation_index] = Array{Tuple{Int32, Float32}}(0)
    for possible_real_tag in keys(probabilities[observation])
      prob_obs_given_real = probabilities[observation][possible_real_tag]
      if (prob_obs_given_real != 0)
        possible_real_tag_index = tag_to_index[possible_real_tag]
        push!(cl[observation_index], (possible_real_tag_index, prob_obs_given_real))
      end
    end
  end
  return cl
end

function tag_index_mapping(tags)
  i = 0
  tag_to_index = Dict{ASCIIString, Int32}()
  index_to_tag = Dict{Int32, ASCIIString}()
  for t in tags
    i += 1
    tag_to_index[t] = i
    index_to_tag[i] = t
  end
  return tag_to_index, index_to_tag
end

function tag_counts(path)
  counts = Dict{ASCIIString, Int32}()
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

function prob_observed_tags_given_reals(counts)
  observed_tags = Set(keys(counts))
  possible_real_tags = observed_tags
  prob_observed_tags_given_reals = Dict{ASCIIString, Dict{ASCIIString, Float32}}()
  for observed_tag in observed_tags
    prob_observed_tags_given_reals[observed_tag] = prob_observed_tag_given_reals(observed_tag, possible_real_tags)
  end
  return prob_observed_tags_given_reals
end

function prob_observed_tag_given_reals(observed_tag, possible_real_tags)
  prob_given_reals = Dict{ASCIIString, Float32}()
  for possible_real_tag in possible_real_tags
    prob_given_reals[possible_real_tag] = 0
  end
  ins_nbrs = insertion_neighbours(observed_tag, possible_real_tags)
  for t in ins_nbrs
    prob_given_reals[t] += ERROR_RATE * INSERTION_RATIO * (1/4)
  end
  del_nbrs = deletion_neighbours(observed_tag, possible_real_tags)
  for t in del_nbrs
    prob_given_reals[t] += ERROR_RATE * DELETION_RATIO
  end
  mut_nbrs = mutation_neighbours(observed_tag, possible_real_tags)
  for t in mut_nbrs
    prob_given_reals[t] += ERROR_RATE * MUTATION_RATIO * (1/3)
  end
  prob_given_reals[observed_tag] = (1 - ERROR_RATE) ^ length(observed_tag)
  return prob_given_reals
end

# function error_neighbours(tag, depth)
#   return error_neighbourhood(tag, depth, 1)
# end
#
# function error_neighbourhood(tag, depth, n)
#   neighbours = Array{ASCIIString}
# end

#The three neighbour functions below can have duplicate tags in the returned lists
#This is used to accumulate probabilities when there are multiple error paths between two tags
# e.g. there are three different deletions that turn AAA into AA

#tag -> insertion -> neighbours
function insertion_neighbours(tag, possible_real_tags)
  neighbours = Array{ASCIIString, 1}()
  word = ""
  for i in 1:length(tag) + 1
    for c in ['A', 'C', 'G', 'T']
      word = insert_at(tag, i, c)
      if (word in possible_real_tags)
        push!(neighbours, word)
      end
    end
  end
  return neighbours
end

#tag -> deletion -> neighbours
function deletion_neighbours(tag, possible_real_tags)
  neighbours = Array{ASCIIString, 1}()
  word = ""
  for i in 1:length(tag)
    word =  without(tag, i)
    if (word in possible_real_tags)
      push!(neighbours, word)
    end
  end
  return neighbours
end

#tag -> mutation -> neighbours
function mutation_neighbours(tag, possible_real_tags)
  neighbours = Array{ASCIIString, 1}()
  word = ""
  for i in 1:length(tag)
    wo = without(tag, i)
    for c in ['A', 'C', 'G', 'T']
      word = insert_at(wo, i, c)
      if (word in possible_real_tags && word != tag)
        push!(neighbours, word)
      end
    end
  end
  return neighbours
end

function without(str, i)
  if i <= 0 || i > length(str)
    return str
  elseif i == 1
    if length(str) > 1
      return str[2:length(str)]
    else
      return ""
    end
  elseif i == length(str)
    return str[1:length(str) - 1]
  else
    return "$(str[1:i-1])$(str[i+1:length(str)])"
  end
end

function insert_at(str, i, c)
  if i <= 1
    return "$c$str"
  elseif i > length(str)
    return "$str$c"
  else
    return "$(str[1:i-1])$c$(str[i:length(str)])"
  end
end

@time process(ARGS[1])
end
