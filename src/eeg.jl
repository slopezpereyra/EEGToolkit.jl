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

    function EEG(file::String; epoch_length::Integer=30, subepoch_length=5, staging::Vector{String}=[""], id::String="")
        id = (id == "") ? file : id

        eeg = EDF.read(file)
        S = Dict()

        for signal in eeg.signals
            x = EDF.decode(signal)
            fs = Int(signal.header.samples_per_record / eeg.header.seconds_per_record)
            S[signal.header.label] = TimeSeries(x, fs; epoch_length=epoch_length, subepoch_length=subepoch_length)
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
`plot_eeg(eeg::EEG, s::Integer, e::Integer; channels::Vector{String}=[""], spacing::AbstractFloat=1.5)`

Plots EEG channels from epoch `s` to epoch `e`. Specific channels may be selected with the `channels` karg.
The `spacing` argument is an added factor in the normalization of the EEG signals - the vertical distance between
each signal in the plot grows proportionally to `spacing`.
"""
function plot_eeg(eeg::EEG, s::Integer, e::Integer; channels::Vector{String}=[""], spacing::AbstractFloat=1.5)

    if channels != [""] && any(x -> x ∉ collect(keys(eeg.signals)), channels)
        throw(ArgumentError("Attempting to plot a non-existent channel. Revise the `channels` keyword argument."))
    end
   
    signal_dict = channels != [""] ? filter(p -> first(p) in channels, eeg.signals) : eeg.signals
    X = [gen_time_domain2(signal.fs, s*signal.epoch_length, (e+1)*signal.epoch_length) for signal in values(signal_dict)]
    L = length(signal_dict)
    p = plot(ylabel = "Amplitude (uV)", xlabel = "Time", yticks = (1:L, collect(keys(signal_dict))));
    for (i, signal) in enumerate(values(signal_dict))
        signal = epoch(signal, s, e).x
        plot!(X[i],  i .+ (signal .- mean(signal)) ./ (spacing * (2*L) * std(signal)), label = "")
    end
    return(p)
end

"""
`remove_channel!(eeg::EEG, channel::String)`

Removes a channel from the EEG.
"""
function remove_channel!(eeg::EEG, channel::String)
    delete!(eeg.signals, channel)
end

"""
`remove_channel!(eeg::EEG, channel::String)`

Removes a list of channels from the EEG.
"""
function remove_channel!(eeg::EEG, channels::Vector{String})
    filter!(p -> first(p) in channels, eeg.signals)
end

"""
`artifact_reject(signal::TimeSeries, anom_matrix::Matrix)`

An anomaly matrix ``A  \\in  \\mathbb{N}^{n  \\times  2}``  is  a  matrix  holds
epoch-subepoch   pairs   that   contain    artifacts    in    a    `TimeSeries`.
Each  row   of   ``(n, m)`` of ``A`` denotes that the ``n``th epoch contained 
an artifact within sub-epoch ``m``. 

This  function  takes  a  `TimeSeries`  signal  and  an  anomaly  matrix  ``A``.
It   removes   from   the   signal   the   epoch-subepoch    pairs    containing
artifacts. In particular, this function returns a 
`Vector{Vector{T}}` with `T<:AbstractFloat`. Thus, if 
`result`   holds   the   return   value   of    this    function,    `result[i]`
contains   what   is   left   from    the    `i`th    epoch    after    removing
its   contaminated   sub-epochs.   It   is   possible   that   `result[i]`    is
empty.
"""
function artifact_reject(signal::TimeSeries, anom_matrix::Matrix)
    T = typeof(signal.x[1])
    epochs = segment(signal, signal.fs * signal.epoch_length)
    windows = map(e -> segment(e, signal.fs * signal.subepoch_length; symmetric=true), epochs)
    for epoch in unique(anom_matrix[:, 1])
        epoch_indexes = findall(x -> x == epoch, anom_matrix[:, 1])
        subeps = anom_matrix[:, 2][epoch_indexes]
        deleteat!(windows[epoch], sort(subeps))
    end
    [Vector{T}( vcat(window...) ) for window in windows]
end


"""
`artifact_reject(signal::TimeSeries, anom_vector::Vector{Integer})`

An anomaly vector ``\\vec{x} \\in \\mathbb{N}^{n}`` is a sorted vector
whose values are those epochs in an `TimeSeries` that contain 
anomalies or artifacts. This function segments the TimeSeries and filters out all 
epochs containing artifacts.
"""
function artifact_reject(signal::TimeSeries, anoms::Vector{Integer})
    T = typeof(signal.x[1])
    epochs = segment(signal, signal.fs * signal.epoch_length)
    deleteat!( epochs, anoms )
end
