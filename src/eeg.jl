using EDF
using Plots
include("ts.jl") 

"""
A struct for the EEG data type.

# Fields
- `signals::Dict{String, TimeSeries}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.
- `staging::Vector{String}`: A vector of stage labels corresponding to each epoch.
- `id::String`: An identifier for the EEG.

# Constructors
`EEG(file::String; epoch_length::Integer=30, staging::Vector{String}=[""], id::String="") -> EEG`: Instantiates an `EEG` from an EDF file (`file`).

# Example
```julia
staging_vector = CSV.read("path/to/stage_data/eeg_staging.csv") # A vector with a stage per each epoch in the EEG
eeg_data = EEG("path/to/edf_data/data.edf"; staging=staging_vector)
```
"""
struct EEG
    signals::Dict{String, TimeSeries}
    staging::Vector{String}
    id::String

    function EEG(file::String; epoch_length::Integer=30, staging::Vector{String}=[""], id::String="")
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
`function filter_by_stage(eeg::EEG, channel::String, stages::Vector) -> `

Returns all portions of an EEG channel in a given stage of the staging vector.
"""
function filter_by_stage(eeg::EEG, channel::String, stages::Vector)
    stage_indexes = findall(x -> x in stages, eeg.staging)
    [epoch(eeg.signals[channel], i) for i in stage_indexes]
end

"""
`artifact_reject(signal::TimeSeries, anom_matrix::Matrix)`

An artifact or anomaly matrix ``A^{n \\times 2}`` is a matrix that marks the time-position of the 
``n`` artifacts in an EEG. Each row of ``A`` is of the form ``(e, s)`` with ``e`` an 
epoch and ``s`` a sub-epoch. This implies an artifact exists in the ``s``th 
sub-epoch of the ``e``th epoch.

This function takes a `TimeSeries` signal and an anomaly matrix ``A``.
It removes from the signal the epoch-subepoch pairs containing 
artifacts. In particular, this function returns a 
`Vector{Vector{T}}` with `T<:AbstractFloat`. Thus, if 
`result` holds the return value of this function, `result[i]`
contains what is left from the `i`th epoch after removing
its contaminated sub-epochs. It is possible that `result[i]` is 
empty.
"""
function artifact_reject(signal::TimeSeries, anom_matrix::Matrix)
    T = typeof(signal.x[1])
    epochs = segment(signal, signal.fs * signal.epoch_length; symmetric=true)
    windows = map(e -> segment(e, signal.fs * 5; symmetric=true), epochs)
    for epoch in unique(anom_matrix[:, 1])
        epoch_indexes = findall(x -> x == epoch, anom_matrix[:, 1])
        subeps = anom_matrix[:, 2][epoch_indexes]
        deleteat!(windows[epoch], sort(subeps))
    end
    [Vector{T}( vcat(window...) ) for window in windows]
end
