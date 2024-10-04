"""
A struct for the EEG data type. An EEG is simply conceived as a collection of labeled time series.

# Fields
- `signals::Dict{String, TimeSeries}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.

# Constructors
`EEG(file::String; id::String="")`: Instantiates an `EEG` from an EDF file (`file`).

# Example
```julia
staging_vector = CSV.read("path/to/stage_data/eeg_staging.csv") # A vector with a stage per each epoch in the EEG
eeg_data = EEG("path/to/edf_data/data.edf"; staging=staging_vector)
```
"""
struct EEG
    _signals::Dict{String, TimeSeries}

    function EEG(file::String)
        if !(endswith(file, ".edf"))
            throw(ArgumentError("The EEG construct must take as input a .edf file, but other filetype was provided."))
        end

        eeg = EDF.read(file)
        S = Dict()

        for signal in eeg.signals
            x = EDF.decode(signal)
            fs = Int(signal.header.samples_per_record / eeg.header.seconds_per_record)
            S[signal.header.label] = TimeSeries(x, fs)
        end

        new(S)
    end
end

"""
example: p -> startswith(first(p), "EEG")
"""
function get_channels(eeg::EEG)::Dict{String, TimeSeries}
    return(eeg._signals)
end


"""
example: p -> startswith(first(p), "EEG")
"""
function get_channel(eeg::EEG, channel_name::String)::TimeSeries
    eeg._signals[channel_name]
end

"""
example: p -> startswith(first(p), "EEG")
"""
function filter_channels(eeg::EEG, filter_function::Function)::Dict{String, TimeSeries}
    filter(filter_function, eeg._signals)
end

"""
example: p -> startswith(first(p), "EEG")
"""
function filter_channels!(eeg::EEG, filter_function::Function)
    filter!(filter_function, eeg._signals)
end


"""
`remove_channel!(eeg::EEG, channel::String)`

Removes a channel from the EEG.
"""
function remove_channel!(eeg::EEG, channel::String)
    delete!(eeg._signals, channel)
end

"""
`remove_channel!(eeg::EEG, channels::Vector{String})`

Removes a list of channels from the EEG.
"""
function remove_channel!(eeg::EEG, channels::Vector{String})
    filter!(p -> first(p) in channels, eeg._signals)
end


"""
`plot_eeg(eeg::EEG, s::Integer, e::Integer; channels::Vector{String}=[""], spacing::AbstractFloat=1.5)`

Plots EEG channels from epoch `s` to epoch `e`. Specific channels may be selected with the `channels` karg.
The `spacing` argument is an added factor in the normalization of the EEG signals - the vertical distance between
each signal in the plot grows proportionally to `spacing`.
"""
function plot_eeg(eeg::EEG, s::Integer, e::Integer; 
                channels::Vector{String}=[""], spacing::AbstractFloat=1.5, 
                epoch_length::Integer=30)

    if channels != [""] && any(x -> x ∉ collect(keys(eeg._signals)), channels)
        throw(ArgumentError("Attempting to plot a non-existent channel. Revise the `channels` keyword argument."))
    end
   
    signal_dict = channels != [""] ? filter(p -> first(p) in channels, eeg._signals) : eeg._signals
    X = [gen_time_domain(signal.fs, s*epoch_length, (e+1)*epoch_length) for signal in values(signal_dict)]
    L = length(signal_dict)
    p = plot(ylabel = "Amplitude (uV)", xlabel = "Time", yticks = (1:L, collect(keys(signal_dict))));
    for (i, signal) in enumerate(values(signal_dict))
        signal = epoch(signal, s, e).x
        plot!(X[i],  i .+ (signal .- mean(signal)) ./ (spacing * (2*L) * std(signal)), label = "")
    end
    return(p)
end
