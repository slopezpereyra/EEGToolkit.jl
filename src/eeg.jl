
"""
A struct for the EEG data type. An EEG is simply conceived as a collection of labeled time series.

# Fields
- `signals::Dict{String, TimeSeries}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.

# Constructors
`EEG(file::String; id::String="")`: Instantiates an `EEG` from an EDF file (`file`).

# Example
```julia
eeg_data = EEG("path/to/edf_data/data.edf")
```
"""
struct EEG
  _signals::Dict{String, TimeSeries}
  _artifacts::Dict{String, ArtifactData}

  function EEG(file::String)
    annotation_flag = false
    if !(endswith(file, ".edf"))
      throw(ArgumentError("The EEG constructor must take as input a .edf file, but other filetype was provided."))
    end

    eeg = EDF.read(file)
    S = Dict()
    A = Dict()

    for signal in eeg.signals
      if signal isa EDF.AnnotationsSignal 
        annotation_flag = true
        continue
      end
      x = EDF.decode(signal)
      fs = Int(signal.header.samples_per_record / eeg.header.seconds_per_record)
      S[signal.header.label] = TimeSeries(x, fs)
      A[signal.header.label] = nothing
    end

    if annotation_flag 
      @warn("The EDF file included annotation signals, which were ignored.")
    end

    new(S, A)
  end
end


"""
get_artifacts(eeg::EEG)::Dict{String, ArtifactData}

Returns artifact data for each channel in the `eeg`.
"""
function get_artifacts(eeg::EEG)::Dict{String, ArtifactData}
  return(eeg._artifacts)
end

"""
get_channels(eeg::EEG)::Dict{String, TimeSeries}

Returns the cannels in the `eeg`.
"""
function get_channels(eeg::EEG)::Dict{String, TimeSeries}
  return(eeg._signals)
end




"""
get_artifacts(eeg::EEG, channel_name::String)::ArtifactData

Returns the artifacts of the channel named `channel_name` from the `eeg`.

"""
function get_artifacts(eeg::EEG, channel_name::String)::ArtifactData
  eeg._artifacts[channel_name]
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
`function detect_artifacts(eeg::EEG, channel_name::String, seg_length::Int)::Nothing`

Performs artifact detection in the `eeg` channel named `channel_name`. 
Stores resulting vector of `Artifact` objects in the `_artifacts` dictionary of the `eeg`
with key `channel_name`. 

Artifact detection is performed on each segment of the signal segmented with
`seg_length`. On each segment, the CAPA algorithm (Fisch et. al 2022) is used 
to detect epidemic distributional changes in the mean value of the segment. 
The `penalty` is an integer value β such that β ln(n) (with `n` the segment's length) 
is the penalty used by the CAPA algorithm to penalize the introduction of artifacts.
The β `penalty` defaults to `24`, which is much higher than the value recommended 
in Fisch et. al but matched human supervision on 78Hz sleep EEGs at the 
developer's laboratory. 
"""
function detect_artifacts(eeg::EEG, channel_name::String, seg_length::Int; penalty::Integer = 24)::Nothing
    
    if isdefined(EEGToolkit, :EEGToolkitR) && isdefined(EEGToolkit.EEGToolkitR, :anomaly)

        signal = get_channel(eeg, channel_name)  
        x = signal.x 
        fs = signal.fs


    start_indices, end_indices, mean_changes, test_statistics = EEGToolkit.EEGToolkitR.anomaly(x, fs, seg_length; penalty)


    artifacts = Artifact[]
    for (start, stop, meanchg, teststat) in zip(start_indices, end_indices, mean_changes, test_statistics)
        epoch = index_to_epoch(start, fs, 30)
        offset = start - (epoch - 1) * (30 * fs)
        subepoch = cld(offset, 5 * fs)

        push!(artifacts, Artifact((start, stop), meanchg, teststat, epoch, subepoch))
    end

        eeg._artifacts[channel_name] = artifacts
        return
    else
        @warn "RCall functionality is not available and artifact detection depends on it. Please install RCall in your current environment."
    end
end


"""
function get_epochs_with_artifacts(eeg::EEG, channel_name::String)

Given an `eeg` and a `channel` that's been artifact detected, returns a vector of all 30-sec epochs which contain an artifact in the channel.
"""
function get_epochs_with_artifacts(eeg::EEG, channel_name::String)
    artifacts = get_artifacts(eeg, channel_name)
    epochs = Int[]
    for a in artifacts
        push!(epochs, a.epoch)  
    end
    return sort(unique(epochs))
end


"""
function plot_artifacts_in_epochs(eeg::EEG, channel_name::String,

Given an `eeg` and an artifact-detected `channel_name`, plots the existing
artifacts in the channel from epoch `from` to epoch `to`. If `annotate` is set
to true, artifacts are annotated with their mean change and test statistic.
"""
function plot_artifacts_in_epochs(eeg::EEG, channel_name::String,
                                  from::Integer, to::Integer; 
                                  annotate::Bool=false)

    signal = get_channel(eeg, channel_name)
    artifacts = get_artifacts(eeg, channel_name)
    plot_artifacts_in_epochs(signal, artifacts, from, to; annotate)
end


