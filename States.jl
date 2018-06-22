push!(LOAD_PATH, ".")

module States
using BioSequences
using Observations

export AbstractState, AbstractStartingState, StartingState
export AbstractRepeatingAnyState, RepeatingAnyState, AbstractObservableState
export ObservableState, AbstractBarcodeState, BarcodeState
export string_to_state_array

abstract type AbstractState end
abstract type AbstractStartingState <: AbstractState end
abstract type AbstractRepeatingAnyState <: AbstractState end
abstract type AbstractObservableState <: AbstractState end
abstract type AbstractBarcodeState <: AbstractObservableState end

struct StartingState <: AbstractStartingState
end

struct RepeatingAnyState <: AbstractRepeatingAnyState
  value::DNA
end

RepeatingAnyState() = RepeatingAnyState(DNA_N)

struct ObservableState <: AbstractObservableState
  value::DNA
end

struct BarcodeState <: AbstractBarcodeState
  value::DNA
end

function string_to_state_array(string_sequence::String)
  states = Array{AbstractState,1}()
  push!(states, StartingState())
  for char in string_sequence
    if char == '*'
      state = RepeatingAnyState()
    elseif islower(char)
      state = BarcodeState(convert(DNA, char))
    else
      state = ObservableState(convert(DNA, char))
    end
    push!(states, state)
  end
  return states
end
end
