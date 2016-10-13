push!(LOAD_PATH, ".")
module Resolving

const REJECT_TAG = "REJECTS"
const ERROR_RATE = 0.005
const DELETION_RATIO = 0.4
const INSERTION_RATIO = 0.4
const MUTATION_RATIO = 0.2

function process(path)
  path = normpath(path)
  counts = tag_counts(path)
  likelihoods = all_likelihood_distributions(counts)
  for k in keys(likelihoods)
    println(k)
    for n in likelihoods[k]
      println(n)
    end
    println("-------------")
  end

  tag_to_index, index_to_tag = tag_index_mapping(Set(keys(counts)))
  println(index_to_tag)

  sparse_probabilities = likelihoods_to_matrix(likelihoods, tag_to_index)
  println(sparse_probabilities)
end

function likelihoods_to_matrix(likelihoods, tag_to_index)
  row_tag_indices = Array{Int64, 1}()
  observed_tag_indices = Array{Int64, 1}()
  probabilities = Array{Float64, 1}()
  all_tags = keys(likelihoods)
  for row_tag in all_tags
    for obs_tag in keys(likelihoods[row_tag])
      p = likelihoods[row_tag][obs_tag]
      if (p != 0)
        rti = tag_to_index[row_tag]
        oti = tag_to_index[obs_tag]
        push!(row_tag_indices, rti)
        push!(observed_tag_indices, oti)
        push!(probabilities, p)
      end
    end
  end
  size = length(all_tags)
  return sparse(row_tag_indices, observed_tag_indices, probabilities, size, size)
end

function tag_index_mapping(tags)
  i = 0
  tag_to_index = Dict{ASCIIString, Int64}()
  index_to_tag = Dict{Int64, ASCIIString}()
  for t in tags
    i += 1
    tag_to_index[t] = i
    index_to_tag[i] = t
  end
  return tag_to_index, index_to_tag
end

function tag_counts(path)
  counts = Dict{ASCIIString, Int64}()
  for f in readdir(path)
    tag, extension = splitext(f)
    if extension == ".fastq" && tag != REJECT_TAG
      lines = 0
      open("$path$f") do tagfile
        lines = countlines(tagfile)
      end
      if (lines % 4 != 0)
        println("Warning! $f appears to have incorrect format: expected file to contain 4 lines per sequence, and to end on a new line.")
      end
      counts[tag] = round(lines / 4)
    end
  end
  return counts
end

function all_likelihood_distributions(counts)
  obs_tags = Set(keys(counts))
  likelihoods = Dict{ASCIIString, Dict{ASCIIString, Float64}}()
  for tag in obs_tags
    likelihoods[tag] = tag_likelihood_distribution(tag, obs_tags)
  end
  return likelihoods
end

function tag_likelihood_distribution(tag, obs_tags)
  likelihood = Dict{ASCIIString, Float64}()
  for obs in obs_tags
    likelihood[obs] = 0
  end
  ins_nbrs = insertion_neighbours(tag)
  for t in ins_nbrs
    if t in obs_tags
      likelihood[t] += ERROR_RATE * INSERTION_RATIO * (1/4)
    end
  end
  del_nbrs = deletion_neighbours(tag)
  for t in del_nbrs
    if t in obs_tags
      likelihood[t] += ERROR_RATE * DELETION_RATIO
    end
  end
  mut_nbrs = mutation_neighbours(tag)
  for t in mut_nbrs
    if t in obs_tags
      likelihood[t] += ERROR_RATE * MUTATION_RATIO * (1/3)
    end
  end
  likelihood[tag] = (1 - ERROR_RATE) ^ length(tag)
  return likelihood
end

#The three neighbour functions below can have duplicate tags in the returned lists
#This is used to accumulate probabilities when there are multiple error paths between two tags
# e.g. there are three different deletions that turn AAA into AA

#tag -> insertion -> neighbours
function insertion_neighbours(tag)
  neighbours = Array{ASCIIString, 1}()
  for i in 1:length(tag) + 1
    for c in ['A', 'C', 'G', 'T']
      push!(neighbours, insert_at(tag, i, c))
    end
  end
  return neighbours
end

#tag -> deletion -> neighbours
function deletion_neighbours(tag)
  neighbours = Array{ASCIIString, 1}()
  for i in 1:length(tag)
    push!(neighbours, without(tag, i))
  end
  return neighbours
end

#tag -> mutation -> neighbours
function mutation_neighbours(tag)
  neighbours = Array{ASCIIString, 1}()
  for i in 1:length(tag)
    wo = without(tag, i)
    for c in ['A', 'C', 'G', 'T']
      rep = insert_at(wo, i, c)
      if (rep != tag)
        push!(neighbours, rep)
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

process(ARGS[1])
end
