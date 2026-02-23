
"""
A struct for the EEG data type. An EEG is conceived as a collection of named
time series with optional masks (named bit vectors) that can be applied to the
time series for various purposes (e.g. artifact removal, sleep staging). Masks
are either global (applicable to all channels) or channel-specific.

# Fields
- `_signals::Dict{String, TimeSeries}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.
- `_global_masks::Dict{Symbol, BitVector}`: A dictionary mapping symbols to bit vectors representing global masks (e.g. NREM masks).
- `_channel_masks::Dict{String, Dict{Symbol, BitVector}}`: A dictionary mapping channel names
to masks (dictionaries mapping symbols to bit vectors). These are channel-specific masks (e.g. artifact masks).

# Constructors
`EEG(file::String; id::String="")`: Instantiates an `EEG` from an EDF file (`file`).

# Example
```julia
eeg_data = EEG("path/to/edf_data/data.edf")
```
"""
struct EEG
    _signals::Dict{String, TimeSeries}
    _global_masks::Dict{Symbol,BitVector}
    _channel_masks::Dict{String,Dict{Symbol,BitVector}}

    function EEG(file::String)

        annotation_flag = false

        endswith(file, ".edf") ||
            throw(ArgumentError("EEG constructor requires .edf file"))

        eeg = EDF.read(file)

        S = Dict{String,TimeSeries}()

        # Global masks start empty
        G = Dict{Symbol,BitVector}()

        # Channel masks
        C = Dict{String,Dict{Symbol,BitVector}}()

        for signal in eeg.signals

            if signal isa EDF.AnnotationsSignal
                annotation_flag = true
                continue
            end

            x = EDF.decode(signal)
            fs = Int(signal.header.samples_per_record /
                     eeg.header.seconds_per_record)

            label = signal.header.label

            S[label] = TimeSeries(x, fs)

            # Initialize empty mask dictionary for this channel
            C[label] = Dict{Symbol,BitVector}()
        end

        if annotation_flag
            @warn("EDF annotation signals were ignored.")
        end

        return new(S, G, C)
    end
end



"""
get_channels(eeg::EEG)::Dict{String, TimeSeries}

Returns the cannels in the `eeg`.
"""
function get_channels(eeg::EEG)::Dict{String, TimeSeries}
  return(eeg._signals)
end

"""
get_channel(eeg::EEG, channel_name::String)::TimeSeries

Returns the channel named `channel_name` from the `eeg`.

"""
function get_channel(eeg::EEG, channel_name::String)::TimeSeries
  eeg._signals[channel_name]
end

"""
filter_channels(eeg::EEG, filter_function::Function)::Dict{String, TimeSeries}

Applies a `filter_function` to the dictionary which maps channel 
names to time series, and returns the filtered result.

Example: `filter_channel(eeg, p -> startswith(first(p), "C"))` would return 
only those EEG channels whose names begin with "C". 

"""
function filter_channels(eeg::EEG, filter_function::Function)::Dict{String, TimeSeries}
  filter(filter_function, eeg._signals)
end

"""
filter_channels!(eeg::EEG, filter_function::Function)::Dict{String, TimeSeries}

Applies by reference a `filter_function` to the dictionary which maps channel 
names to time series. 

Example: `filter_channel(eeg, p -> startswith(first(p), "C"))` would keep 
only those EEG channels whose names begin with "C". 

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


"""
    add_mask!(eeg::EEG, name::Symbol, mask::BitVector)

Adds a mask with key `name` and value `mask` to the global masks of `eeg`. Global masks apply to all channels (e.g. a NREM masks).
"""
function add_mask!(eeg::EEG, name::Symbol, mask::BitVector)
    eeg._global_masks[name] = mask
end


"""
    add_mask!(eeg::EEG, channel::String, name::Symbol, mask::BitVector)

Adds a mask with key `name` and value `mask` to the masks of `channel` in `eeg`.
Channel-specific masks apply only to a single channel (e.g. an artifact mask).
"""

function add_mask!(eeg::EEG, channel::String, name::Symbol, mask::BitVector)
    eeg._channel_masks[channel][name] = mask
end


"""
    get_masks(eeg::EEG)::Dict{Symbol,BitVector}

Return all global (EEG-level) masks associated with `eeg`.

Global masks apply to all channels (e.g. sleep staging, NREM masks).
"""
function get_masks(eeg::EEG)::Dict{Symbol,BitVector}
    return eeg._global_masks
end


"""
    get_masks(eeg::EEG, channel::String)::Dict{Symbol,BitVector}

Return all masks associated with `channel`.

This includes only channel-specific masks (e.g. artifact masks).
"""
function get_masks(eeg::EEG, channel::String)::Dict{Symbol,BitVector}

    haskey(eeg._channel_masks, channel) ||
        throw(ArgumentError("Channel '$channel' not found."))

    return eeg._channel_masks[channel]
end


"""
    get_masks(eeg::EEG, name::Symbol)::BitVector

Return a global mask by name.
"""
function get_masks(eeg::EEG, name::Symbol)::BitVector

    haskey(eeg._global_masks, name) ||
        throw(ArgumentError("Global mask :$name not found."))

    return eeg._global_masks[name]
end


"""
    get_masks(eeg::EEG, channel::String, name::Symbol)::BitVector

Return mask `name` for the given channel.
"""
function get_masks(eeg::EEG, channel::String, name::Symbol)::BitVector

    haskey(eeg._channel_masks, channel) ||
        throw(ArgumentError("Channel '$channel' not found."))

    masks = eeg._channel_masks[channel]

    haskey(masks, name) ||
        throw(ArgumentError("Mask :$name not found for channel '$channel'."))

    return masks[name]
end


"""
    describe(eeg::EEG; epoch_length::Int=30)

Prints a description of the `EEG` object and returns a `DataFrame` with details 
specific to each channel. 
"""
function describe(eeg::EEG; epoch_length::Int=30)
    # Extraemos las señales utilizando tu API
    signals = get_channels(eeg)
    n_channels = length(signals)
    
    println("="^45)
    println("EEG Summary")
    println("="^45)
    println("Total Channels: ", n_channels)
    
    global_masks = get_masks(eeg)
    println("Global Masks:   ", length(global_masks))
    if !isempty(global_masks)
        println("  ↳ ", join(collect(keys(global_masks)), ", "))
    end
    println("-"^45)

    if n_channels == 0
        return DataFrame()
    end
    
    # Sintaxis segura para nombres de columnas con espacios/caracteres especiales
    df = DataFrame(
        "Channel" => String[], 
        "Fs (Hz)" => Int[], 
        "Duration" => Time[], 
        "Epochs" => Int[],
        "Masks" => Int[]
    )
    
    for (name, ts) in signals
        sec = length_in_secs(ts)
        ep = num_epochs(ts, epoch_length)
        c_masks = length(get_masks(eeg, name)) 
        
        # Aprovechamos tu función seconds_to_time de src/ts.jl
        dur_time = seconds_to_time(sec)
        
        push!(df, (name, ts.fs, dur_time, ep, c_masks))
    end
    
    return df
end
