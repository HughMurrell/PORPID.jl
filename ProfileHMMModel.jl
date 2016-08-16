push!(LOAD_PATH, ".")
module ProfileHMMModel

using Nucleotides
using Observations

export string_to_state_model, viterbi
end
