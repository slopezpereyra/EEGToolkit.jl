using EDF
using Plots


mutable struct EEG
    """
    A mutable struct representing EEG data with associated metadata.

    # Fields
    - `signals::Dict{String, Vector{<:AbstractFloat}}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.
    - `sampling_rates::Dict{String, Integer}`: A dictionary mapping signals (strings) to the integer sampling rates.
    - `fs::Integer`: Sampling frequency (in Hz) of the EEG signals.
    - `epoch_length::Integer`: Length of each epoch (in seconds).
    - `staging::Vector{String}`: A vector of stage labels corresponding to each epoch.
    - `id::String`: An optional identifier for the EEG.

    # Constructor
    - `EEG(file, fs, epoch_length, staging)`: Constructs an `EEG` object from an EDF file (`file`) containing EEG data. The function reads the signals, computes the necessary metadata (`fs`, `N`, `epoch_count`), and initializes the `EEG` struct with the provided `staging` vector.

    ## Example
    ```julia
    eeg_data = EEG("data.edf", 256, 30, ["Wake", "1", "1", …, "REM", "REM", "2", "3", …])
    """
    signals::Dict{String,Vector{<:AbstractFloat}}
    sampling_rates::Dict{String,Integer}
    fs::Integer
    epoch_length::Integer
    staging::Vector{String}
    id::String

    function EEG(file, epoch_length, staging, id="")

        eeg = EDF.read(file)
        S = Dict()
        FS = Dict()

        for signal in eeg.signals
            S[signal.header.label] = EDF.decode(signal)
            FS[signal.header.label] = Int(signal.header.samples_per_record / eeg.header.seconds_per_record)
        end

        new(S, FS, first(values(FS)), epoch_length, staging, id)
    end

end


"""
An epoch e is a function ℕ ↦ ℕ₀² s.t. e(n) = (x, y) if and only if 
[x, y] is the space of all index values in the ith epoch of an EEG signal.

Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
[x, y] is the space of all indexes corresponding to values in an EEG signal 
in epochs n, n +1, …, m.
 """
function epoch(eeg::EEG, n::Integer)
    return [((n - 1) * eeg.fs * eeg.epoch_length) + 1, n * eeg.fs * eeg.epoch_length]
end

function epoch(eeg::EEG, n::Integer, m::Integer)
    """An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
       [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

        Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
        [x, y] is the space of all indexes corresponding to values in an EEG signal 
        in epochs n, n +1, …, m.
       """
    if (n == m)
        return epoch(eeg, n)
    end
    if (n > m)
        throw(ArgumentError("The second epoch should be greater than the first."))
    end
    return [((n - 1) * eeg.fs * eeg.epoch_length) + 1, m * eeg.fs * eeg.epoch_length]
end

"""
    An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
    [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

    Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
    [x, y] is the space of all indexes corresponding to values in an EEG signal 
    in epochs n, n +1, …, m.
"""
function epoch(eeg::EEG, n::Integer, channel::String)
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n)
    signal[bounds[1]:bounds[2]]
end

"""
    An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
    [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

    Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
    [x, y] is the space of all indexes corresponding to values in an EEG signal 
    in epochs n, n +1, …, m.
"""
function epoch(eeg::EEG, n::Integer, m::Integer, channel::String)
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n, m)
    signal[bounds[1]:bounds[2]]
end

"""Generates time domain of an EEG signal from second `s` to second `e`."""
function time(eeg::EEG, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})
    step = 1 / eeg.fs
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
    t = time(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signal = epoch(eeg, n, m, channels)
    plot(t, signal)
end

"""
    Plots with the active backend all EEG channels in `channels` in the 
    range from the `n`th to the `m`th epoch.
"""
function plot_eeg(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    plot(t, signals, layout=(length(signals), 1), legend=false)
    xlabel!("Time")
    ylabel!("")
end

"""
    Plots with the active backend all EEG channels in `channels` in the 
    range from the `n`th to the `m`th epoch.
"""
function plot_eeg_overlay(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    Plots.plot(t, signals, legend=false)
    xlabel!("Time")
    ylabel!("")
end


"""
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
   Returns the subset of an EEG channel corresponding to a stage. 
"""
function get_stage(eeg::EEG, channel::String, stages::Vector)
    indexes = get_stage_indexes(eeg, stages)
    eeg.signals[channel][indexes]
end

function sfilter!(eeg::EEG, channel::String, digfilter, cut_off)
    eeg.signals[channel] = filt(digitalfilter(digfilter(cut_off, fs=eeg.fs), Butterworth(4)), eeg.signals[channel])
end

function sfilter!(eeg::EEG, channels::Vector{<:String}, digfilter, cut_off)

    for chan in channels
        eeg.signals[chan] = filt(digitalfilter(digfilter(cut_off, eeg.fs), Butterworth(4)), eeg.signals[chan])
    end
end

function sfilter!(eeg::EEG, digfilter, cut_off)
    for chan in keys(eeg.signals)
        eeg.signals[chan] = filt(digitalfilter(digfilter(cut_off, eeg.fs), Butterworth(4)), eeg.signals[chan])
    end
end

function artifact_reject(eeg, anom_matrix, signal)
    x = eeg.signals[signal]
    epochs = overlaps(x, eeg.fs*30, 0)
    windows = map(x -> overlaps(x, eeg.fs * 5, 1/(eeg.fs*5)), epochs)
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










