push!(LOAD_PATH, ".")

module States
using Nucleotides
using Observations
export ObservableState, StartingState, RepeatingAnyState, State, string_to_state_model

type ObservableState
  value::DNASymbol
end

type StartingState
end

type RepeatingAnyState
end

State = Union{ObservableState, StartingState, RepeatingAnyState}

function string_to_state_model(string_sequence::AbstractString)
  states = Array{State,1}()
  push!(states, StartingState())
  for char in string_sequence
    if char == '*'
      state = RepeatingAnyState()
    else
      state = ObservableState(convert(DNASymbol, char))
    end
    push!(states, state)
  end
  return states
end
end
