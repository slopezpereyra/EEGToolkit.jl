using EDF
using Plots


"""
A mutable struct representing EEG data with associated metadata.

# Fields
- `signals::Dict{String, Vector{<:AbstractFloat}}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.
- `sampling_rates::Dict{String, Integer}`: A dictionary mapping signals (strings) to the integer sampling rates.
- `fs::Integer`: A default sampling rate that will be used for calculations. Defaults to the maximum sampling rate among all signals.
- `epoch_length::Integer`: Length of each epoch (in seconds).
- `staging::Vector{String}`: A vector of stage labels corresponding to each epoch.
- `id::String`: An identifier for the EEG.

# Constructors
`EEG(file, fs, epoch_length, staging)`: Constructs an `EEG` object from an EDF file (`file`) containing EEG data. The function reads the signals, computes the necessary metadata (`fs`, `N`, `epoch_count`), and initializes the `EEG` struct with the provided `staging` vector.

# Example
```julia
staging_vector = CSV.read("path/to/stage_data/eeg_staging.csv") # A vector with a stage per each epoch in the EEG
eeg_data = EEG("path/to/edf_data/data.edf", 30, staging_vector)

# Alternatively, if no stage data exists, it is safe to do 
eeg_data = EEG("path/to/edf_data/data.edf", 30, [])

# or simply 
eeg_data = EEG("path/to/edf_data/data.edf", 30) # Staging vector defaults to [] if not provided

# or simply
eeg_data = EEG("path/to/edf_data/data.edf") # Epoch length defaults to 30 if not provided.
```
"""
mutable struct EEG
    signals::Dict{String,Vector{<:AbstractFloat}}
    sampling_rates::Dict{String,Integer}
    fs::Integer
    epoch_length::Integer
    staging::Vector{String}
    id::String

    function EEG(file::String, epoch_length::Integer=30, staging::Vector{String}=[""], id::String="")
        id = (id == "") ? file : id

        eeg = EDF.read(file)
        S = Dict()
        FS = Dict()

        for signal in eeg.signals
            S[signal.header.label] = EDF.decode(signal)
            FS[signal.header.label] = Int(signal.header.samples_per_record / eeg.header.seconds_per_record)
        end

        fₛ = maximum(values(FS))

        new(S, FS, fₛ, epoch_length, staging, id)
    end
end


"""
`epoch(eeg::EEG, n::Integer, fs::Integer=-1)`

Returns a vector [i₁, …, iₖ] with all indexes corresponding to the `n`th epoch of the EEG.
The default sampling rate is used to compute the indexes.
"""
function epoch(eeg::EEG, n::Integer, fs::Integer=-1)
    fs = (fs == -1) ? eeg.fs : fs
    return [((n - 1) * fs * eeg.epoch_length) + 1, n * fs * eeg.epoch_length]
end

"""
`epoch(eeg::EEG, n::Integer, m::Integer, fs::Integer=-1)`

Returns a vector [i₁, …, iₖ] with all indexes corresponding to epochs `n, n+1, …, m` of the EEG.
The default sampling rate is used to compute the indexes.
"""
function epoch(eeg::EEG, n::Integer, m::Integer, fs::Integer=-1)
    if (n == m)
        return epoch(eeg, n)
    end
    if (n > m)
        throw(ArgumentError("The second epoch should be greater than the first."))
    end
    fs = (fs == -1) ? eeg.fs : fs
    return [((n - 1) * eeg.fs * eeg.epoch_length) + 1, m * eeg.fs * eeg.epoch_length]
end

"""
`epoch(eeg::EEG, n::Integer, channel::String)`

Returns a vector [x₁, …, xₖ] with all values of the signal `channel` in the `n`th epoch.
"""
function epoch(eeg::EEG, n::Integer, channel::String)
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n, eeg.sampling_rates[channel])
    signal[bounds[1]:bounds[2]]
end

"""
`epoch(eeg::EEG, n::Integer, m::Integer, channel::String)`

Returns a vector [x₁, …, xₖ] with all values of the signal `channel` in the epochs `n, n+1, …, m`.
"""
function epoch(eeg::EEG, n::Integer, m::Integer, channel::String)
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n, m, eeg.sampling_rates[channel])
    signal[bounds[1]:bounds[2]]
end

"""
`gen_time_domain(eeg::EEG, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer}, fs::Integer=-1)`

Given an EEG, generates the time vector t₁, …, tₙ corresponding to 
EEG signals from time `s` to `e`.
"""
function gen_time_domain(eeg::EEG, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer}, fs::Integer=-1)
    fs = (fs == -1) ? eeg.fs : fs
    step = 1 / fs
    return [i for i in (s+step):step:e]
end

"""
Plots with the active backend all EEG channels in `channels` in the 
range from the `n`th to the `m`th epoch.
"""
function plot_eeg(eeg::EEG, channels::String, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = gen_time_domain(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signal = epoch(eeg, n, m, channels)
    plot(t, signal)
end

"""
`plot_eeg(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)`

Plots with the active backend all EEG channels in `channels` in the 
range from the `n`th to the `m`th epoch.
"""
function plot_eeg(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = gen_time_domain(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    plot(t, signals, layout=(length(signals), 1), legend=false)
    xlabel!("Time")
    ylabel!("")
end

"""
`plot_eeg_overlay(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)`

Plots with the active backend all EEG channels in `channels` in the 
range from the `n`th to the `m`th epoch.
"""
function plot_eeg_overlay(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = gen_time_domain(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    Plots.plot(t, signals, legend=false)
    xlabel!("Time")
    ylabel!("")
end


"""
`get_stage_indexes(eeg::EEG, stages::Vector)`

This function maps an EEG and a`stages` vector to the array of all indexes whose values in an EEG 
signal pertain to a stage in `stages`.
"""
function get_stage_indexes(eeg::EEG, stages::Vector)

    stage_indexes = findall(x -> x in stages, eeg.staging)
    result = []

    for index in stage_indexes
        epoch_range = epoch(eeg, index)
        v = collect(epoch_range[1]:epoch_range[2])
        result = vcat(result, v)
    end

    return result
end

"""
`function get_stage(eeg::EEG, channel::String, stages::Vector)`

Returns all portions of an EEG channel in a given stage of the staging vector.
"""
function get_stage(eeg::EEG, channel::String, stages::Vector)
    indexes = get_stage_indexes(eeg, stages)
    eeg.signals[channel][indexes]
end

"""
`artifact_reject(eeg::EEG, anom_matrix::Matrix, signal::String)`

Given an EEG, a 2x2 matrix associating epoch-subepoch pairs with artifacts, and a signal,
returns a subset of the signal with all sub-epochs containing artifacts removed.

The signal is split in epoch-length windows and each window is split in subepoch-length 
windows; the matrix gives the epoch and subepoch indexes to be removed. 
"""
function artifact_reject(eeg::EEG, anom_matrix::Matrix, signal::String)
    x = eeg.signals[signal]
    epochs = overlaps(x, eeg.fs * 30, 0)
    windows = map(x -> overlaps(x, eeg.fs * 5, 1 / (eeg.fs * 5)), epochs)
    for epoch in unique(anom_matrix[:, 1])
        if epoch > length(windows)
            @warn("Anomaly was found on final epoch, which does not belong to the segmentation\n")
            continue
        end
        subeps_indexes = findall(x -> x == epoch, anom_matrix[:, 1])
        subeps = anom_matrix[:, 2][subeps_indexes]
        deleteat!(windows[epoch], sort(subeps))
    end
    clean = map(x -> vcat(x...), windows)
    return clean
end

