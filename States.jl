push!(LOAD_PATH, ".")

module States
using Nucleotides
using Observations

export AbstractState, AbstractStartingState, StartingState
export AbstractRepeatingAnyState, RepeatingAnyState, AbstractObservableState
export ObservableState, AbstractBarcodeState, BarcodeState
export string_to_state_array

abstract AbstractState
abstract AbstractStartingState <: AbstractState
abstract AbstractRepeatingAnyState <: AbstractState
abstract AbstractObservableState <: AbstractState
abstract AbstractBarcodeState <: AbstractObservableState

type StartingState <: AbstractStartingState
end

type RepeatingAnyState <: AbstractRepeatingAnyState
end

type ObservableState <: AbstractObservableState
  value::DNASymbol
end

type BarcodeState <: AbstractBarcodeState
  value::DNASymbol
end

function string_to_state_array(string_sequence::AbstractString)
  states = Array{AbstractState,1}()
  push!(states, StartingState())
  for char in string_sequence
    if char == '*'
      state = RepeatingAnyState()
    elseif islower(char)
      state = BarcodeState(convert(DNASymbol, char))
    else
      state = ObservableState(convert(DNASymbol, char))
    end
    push!(states, state)
  end
  return states
end
end
