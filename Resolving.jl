push!(LOAD_PATH, ".")
module Resolving

using SparseICMapLDA

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
  @time likelihoods = all_likelihood_distributions(counts)

  println(STDERR, "Generating index mapping...")
  @time tag_to_index, index_to_tag = tag_index_mapping(Set(keys(counts)))

  indexed_counts = index_counts(counts, tag_to_index)

  println(STDERR, "Converting tag mapping to sparse matrix...")
  @time sparse_probabilities = likelihoods_to_matrix(likelihoods, tag_to_index)

  println(STDERR, "Iterating...")
  @time posterior = SparseICMapLDA.LDA(sparse_probabilities, indexed_counts)

  println("tag,count,best_tag,best_score")
  n = size(posterior, 1)
  for r in 1:n
    maxval = -1
    maxindex = 0
    for c in 1:n
      if posterior[r,c] > maxval
        maxval = posterior[r,c]
        maxindex = c
      end
    end
    #println("$(index_to_tag[r])\t$(round(posterior[r,maxindex], 4))")
    #if (maxval < 0.99)
    println("$(index_to_tag[r]),$(counts[index_to_tag[r]]),$(index_to_tag[maxindex]),$(round(posterior[r,maxindex], 4))")
    #end
  end
end

function index_counts(counts, tag_to_index)
  tags = keys(counts)
  ic = zeros(Int32, length(tags), 1)
  for t in tags
    ic[tag_to_index[t], 1] = counts[t]
  end
  return ic
end

function likelihoods_to_matrix(likelihoods, tag_to_index)
  row_tag_indices = Array{Int32, 1}()
  observed_tag_indices = Array{Int32, 1}()
  probabilities = Array{Float32, 1}()
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

function all_likelihood_distributions(counts)
  obs_tags = Set(keys(counts))
  likelihoods = Dict{ASCIIString, Dict{ASCIIString, Float32}}()
  for tag in obs_tags
    likelihoods[tag] = tag_likelihood_distribution(tag, obs_tags)
  end
  return likelihoods
end

function tag_likelihood_distribution(tag, obs_tags)
  likelihood = Dict{ASCIIString, Float32}()
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

@time process(ARGS[1])
end
