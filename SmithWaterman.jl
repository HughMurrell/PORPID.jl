push!(LOAD_PATH, ".")
module SmithWaterman
using Nucleotides
using Observations
using States
export extract_tag

const STAR_INSERTION_SCORE = 0

@enum AlignOp OP_MATCH=1 OP_DEL=2 OP_INS=3

function extract_tag(observations::Array{Observation,1}, states::Array{State,1})
  rows = length(states)
  cols = length(observations) + 1 #first column is for 'before first symbol' position
  scores = zeros(Float64, rows, cols)
  ops = Array{AlignOp}(rows, cols)
  #Left-most column contains all deletions
  for r = 2:rows
    if typeof(states[r]) <: RepeatingAnyState
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
      if typeof(currentstate) <: RepeatingAnyState
        #RepeatingAnyState represents any number of insertions at star_insertion_score or can be freely deleted
        scores[r,c] = max(scores[r-1,c], scores[r,c-1] + STAR_INSERTION_SCORE)
      else
        matchscore = scores[r-1,c-1] + log(prob(currentstate.value, observations[c-1].value, observations[c-1].prob))
        inscore = scores[r,c-1] + L_PROBABILITY_OF_INSERTION
        delscore = scores[r-1,c] + L_PROBABILITY_OF_DELETION
        #The order of operations is significant in the case of multiple paths with the same score
        #In this case, a single path is chosen based on operation order
        if (delscore >= matchscore && delscore >= inscore)
          scores[r,c] = delscore
          ops[r,c] = OP_DEL
        elseif (inscore >= matchscore)
          scores[r,c] = inscore
          ops[r,c] = OP_INS
        else
          scores[r,c] = matchscore
          ops[r,c] = OP_MATCH
        end
      end
    end
  end
  #Construct tag from score and operation matrices
  tag = Array{DNASymbol, 1}()
  r = rows
  c = cols
  while (r > 1 || c > 1)
    if (ops[r,c] == OP_MATCH)
      if (states[r].value == Nucleotides.DNA_N)
        insert!(tag, 1, observations[c-1].value)
      end
      r = r - 1
      c = c - 1
    elseif (ops[r,c] == OP_INS)
      c = c - 1
    else
      r = r - 1
    end
  end
  return (scores[rows,cols], tag)
end
end
