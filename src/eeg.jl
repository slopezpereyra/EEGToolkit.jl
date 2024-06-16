using EDF
using Plots
include("ts.jl") 

"""
A struct representing EEG data with associated metadata.

# Fields
- `signals::Dict{String, TimeSeries}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.
- `staging::Vector{String}`: A vector of stage labels corresponding to each epoch.
- `id::String`: An identifier for the EEG.

# Constructors
`EEG(file::String, epoch_length::Integer=30, staging::Vector{String}=[""], id::String="")`: Constructs an `EEG` object from an EDF file (`file`) containing EEG data. The function reads the signals, computes the necessary metadata (`fs`, `N`, `epoch_count`), and initializes the `EEG` struct with the provided `staging` vector.

# Example
```julia
staging_vector = CSV.read("path/to/stage_data/eeg_staging.csv") # A vector with a stage per each epoch in the EEG
eeg_data = EEG("path/to/edf_data/data.edf", 30, staging_vector)

# Alternatively, if no stage data exists, it is safe to do 
eeg_data = EEG("path/to/edf_data/data.edf", 30) # Staging vector defaults to [] if not provided

# or simply
eeg_data = EEG("path/to/edf_data/data.edf") # Epoch length defaults to 30 if not provided.
```
"""
struct EEG
    signals::Dict{String, TimeSeries}
    staging::Vector{String}
    id::String

    function EEG(file::String, epoch_length::Integer=30, staging::Vector{String}=[""], id::String="")
        id = (id == "") ? file : id

        eeg = EDF.read(file)
        S = Dict()

        for signal in eeg.signals
            x = EDF.decode(signal)
            fs = Int(signal.header.samples_per_record / eeg.header.seconds_per_record)
            S[signal.header.label] = TimeSeries(x, fs, epoch_length)
        end

        new(S, staging, id)
    end
end


"""
`function filter_by_stage(eeg::EEG, channel::String, stages::Vector)`

Returns all portions of an EEG channel in a given stage of the staging vector.
"""
function filter_by_stage(eeg::EEG, channel::String, stages::Vector)
    stage_indexes = findall(x -> x in stages, eeg.staging)
    [epoch(eeg.signals[channel], i) for i in stage_indexes]
end

"""
`artifact_reject(signal::TimeSeries, anom_matrix::Matrix)`

Given an TimeSeries and a 2x2 matrix associating epoch-subepoch pairs with artifacts, 
returns a new TimeSeries with all sub-epochs containing artifacts removed.
"""
function artifact_reject(signal::TimeSeries, anom_matrix::Matrix)
    epochs = segment(signal, signal.fs * signal.epoch_length; symmetric=true)
    windows = map(e -> segment(e, signal.fs * 5; symmetric=true), epochs)
    for epoch in unique(anom_matrix[:, 1])
        epoch_indexes = findall(x -> x == epoch, anom_matrix[:, 1])
        subeps = anom_matrix[:, 2][epoch_indexes]
        deleteat!(windows[epoch], sort(subeps))
    end
    [vcat(window...) for window in windows]
end
