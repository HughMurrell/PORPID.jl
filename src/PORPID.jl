module PORPID

export # States
    AbstractState,
    AbstractStartingState,
    StartingState,
    AbstractRepeatingAnyState,
    RepeatingAnyState,
    AbstractObservableState,
    ObservableState,
    AbstractBarcodeState,
    BarcodeState,
    string_to_state_array
export # LDA
    LDA
export # NeedlemanWunsch
    extract_tag
export # Observations
    Observation,
    phred_score_to_prob,
    prob,
    PROBABILITY_OF_INSERTION,
    PROBABILITY_OF_DELETION,
    L_PROBABILITY_OF_INSERTION,
    L_PROBABILITY_OF_DELETION,
    L_PROBABILITY_PER_EXTRA_BASE,
    L_PROBABILITY_OF_NORMAL_TRANSITION
export # PORPIDConfig
    Configuration,
    Template,
    fasta,
    fastq,
    load_config_from_json,
    extract_tags,
    extract_tags_from_file,
    sequence_to_observation,
    best_of_forward_and_reverse,
    slice_sequence
export # Resolvng
    resolve_tags_in_dir,
    tag_index_mapping,
    prob_observed_tags_given_reals,
    index_counts

using BioSequences

include("States.jl")
include("CustomLDA.jl")
include("NeedlemanWunsch.jl")
include("Observations.jl")
include("PORPIDConfig.jl")
include("PORPIDMethods.jl")
include("Resolving.jl")

end # module
