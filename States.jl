push!(LOAD_PATH, ".")

module States
using Bio.Seq
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
  value::DNANucleotide
end

RepeatingAnyState() = RepeatingAnyState(DNA_N)

type ObservableState <: AbstractObservableState
  value::DNANucleotide
end

type BarcodeState <: AbstractBarcodeState
  value::DNANucleotide
end

function string_to_state_array(string_sequence::String)
  states = Array{AbstractState,1}()
  push!(states, StartingState())
  for char in string_sequence
    if char == '*'
      state = RepeatingAnyState()
    elseif islower(char)
      state = BarcodeState(convert(DNANucleotide, char))
    else
      state = ObservableState(convert(DNANucleotide, char))
    end
    push!(states, state)
  end
  return states
end
end
