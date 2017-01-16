push!(LOAD_PATH, ".")
module NeedlemanWunsch
using Nucleotides
using Observations
using States
export extract_tag

const STAR_INSERTION_SCORE = 0
const SCORE_EPSILON = 1e-10
#SCORE_EPSILON is used to prevent accuracy limitations from inducing false inequalities
#e.g. log(0.25) + log(0.01) + log(0.25) != log(0.25) + log(0.25) + log(0.01)

@enum AlignOp OP_MATCH=1 OP_DEL=2 OP_INS=3

function extract_tag(observations::Array{Observation,1}, states::Array{AbstractState,1})
  rows = length(states)
  cols = length(observations) + 1 #first column is for 'before first symbol' position
  scores = zeros(Float64, rows, cols)
  ops = Array{AlignOp}(rows, cols)
  #Left-most column contains all deletions
  for r = 2:rows
    if typeof(states[r]) <: AbstractRepeatingAnyState
      scores[r,1] = scores[r-1,1] #Deletion of RepeatingAnyState is free
    else
      scores[r,1] = scores[r-1,1] + L_PROBABILITY_OF_DELETION
    end
    ops[r,1] = OP_DEL
  end
  #Top row contains all insertions
  for c = 2:cols
    scores[1,c] = scores[1,c-1] + L_PROBABILITY_OF_INSERTION
    ops[1,c] = OP_INS
  end
  #Fill in scores for the rest of the matrix
  for r = 2:rows
    for c = 2:cols
      currentstate = states[r]
      if typeof(currentstate) <: AbstractRepeatingAnyState
        #RepeatingAnyState represents any number of insertions at star_insertion_score or can be freely deleted
        inscore = scores[r,c-1] + STAR_INSERTION_SCORE
        delscore = scores[r-1,c]
        if delscore > inscore - SCORE_EPSILON
          scores[r,c] = delscore
          ops[r,c] = OP_DEL
        else
          scores[r,c] = inscore
          ops[r,c] = OP_INS
        end
      else
        matchscore = scores[r-1,c-1] + log(prob(currentstate.value, observations[c-1].value, observations[c-1].prob))
        inscore = scores[r,c-1] + L_PROBABILITY_OF_INSERTION
        delscore = scores[r-1,c] + L_PROBABILITY_OF_DELETION
        #The order of operations is significant in the case of multiple paths with the same score
        #In this case, a single path is chosen based on operation order
        #Currently matching is done first (when backtracing) so that insertions in the tag section are aligned to N symbols
        #If early insertion is preferred over early matching, such insertions might be aligned to known symbols in the template instead
        if matchscore > delscore - SCORE_EPSILON && matchscore > inscore - SCORE_EPSILON
          scores[r,c] = matchscore
          ops[r,c] = OP_MATCH
        elseif inscore > delscore - SCORE_EPSILON
          scores[r,c] = inscore
          ops[r,c] = OP_INS
        else
          scores[r,c] = delscore
          ops[r,c] = OP_DEL
        end
      end
    end
  end

  #Construct tag from score and operation matrices
  tag = Array{DNASymbol, 1}()
  errors = 0
  r = rows
  c = cols
  while r > 1 || c > 1
    #Debugging print: path and operations
    #print("r=$r\tc=$c\t$(ops[r,c])$(ops[r,c]!=OP_MATCH ? "\t\t" : "\t")Obs=$(c > 1 ? string(observations[c-1].value) : "0")")
    #print("\tState=$(typeof(states[r]) <: AbstractRepeatingAnyState ? "*" : (typeof(states[r]) <: AbstractStartingState ? "0" : string(states[r].value)))")
    #print("\t\tScore=$(round(scores[r,c], 2))\n")
    if ops[r,c] == OP_MATCH
      if typeof(states[r]) <: AbstractBarcodeState
        insert!(tag, 1, observations[c-1].value)
      end
      if typeof(states[r]) <: AbstractObservableState
        if !(observations[c-1].value in Nucleotides.constituents(states[r].value))
          errors = errors + 1
        end
      end
      r = r - 1
      c = c - 1
    elseif ops[r,c] == OP_INS
      #Include in the tag any insertions that are aligned to an N or aligned adjacent to an N
      #No need to check states[r-1] as long as matching is done before insertion, because any insertions to the right of a series of Ns will
      #  be matched to the Ns and the symbols aligned to the left of the Ns will be considered insertions instead (and included in the tag)
      if typeof(states[r]) <: AbstractBarcodeState ||
          (r < rows && typeof(states[r+1]) <: AbstractBarcodeState)
        insert!(tag, 1, observations[c-1].value)
      end
      if !(typeof(states[r]) <: AbstractRepeatingAnyState)
        errors = errors + 1
      end
      c = c - 1
    else #OP_DEL
      if !(typeof(states[r]) <: AbstractRepeatingAnyState)
        errors = errors + 1
      end
      r = r - 1
    end
  end
  return (scores[rows,cols], tag, errors)
end
end
