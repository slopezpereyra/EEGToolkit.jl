using CSV
using DataFrames
using Plots
using RCall

R"""
library(dplyr)
library(anomaly)
library(ggplot2)
library(purrr)
library(tidyr)

artifact_detection <- function(signal, fs, epoch_length_in_secs){
    epoch_length <- epoch_length_in_secs * fs
    n_epochs <- floor(length(signal) / epoch_length)
    epochs <- split(signal[1:(n_epochs * epoch_length)],
                    rep(1:n_epochs, each = epoch_length))
    f <- function(x) { anomaly::capa(scale(x), type="mean", beta=24*log(length(x))) }
    analyses <- purrr::map(epochs, f)
    return(analyses)
}

get_global_artifacts <- function(i, analyses, epoch_length, fs) {
  offset <- (i - 1) * epoch_length * fs
  canoms <- anomaly::collective_anomalies(analyses[[i]])
  if (nrow(canoms) == 0) return(NULL)
  canoms$start <- canoms$start + offset
  canoms$end   <- canoms$end + offset
  canoms$epoch <- i
  return(canoms)
}

detect_artifacts <- function(signal, fs, epoch_length) {
    analyses <- artifact_detection(signal, fs, epoch_length)
    g <- function(i) { get_global_artifacts(i, analyses, epoch_length, fs) }
    global_artifacts <- as_tibble(purrr::map_dfr(seq_along(analyses), g))
    global_artifacts <- global_artifacts %>%
        filter(test.statistic > 1000) %>%
        select(start, end, mean.change, test.statistic, epoch) %>%
        mutate(start = as.integer(start),
               end = as.integer(end))
    return(global_artifacts)
}
"""

struct Artifact
  loc::Tuple{Int, Int}
  mean_change::Float64 
  test_statistic::Float64
  epoch_range::Tuple{Int, Int}
end


"""
detect_artifacts(signal::TimeSeries, epoch_length::Int)::DataFrame

Detects collective artifacts in a time series using the CAPA algorithm from 
the R `artifact` package. The signal is divided into non-overlapping epochs 
of `epoch_length` seconds, each of which is independently analyzed for 
mean-shift artifacts. The results are combined into a single `DataFrame` 
containing the global start and end sample indices of each artifact, 
the mean change, the test statistic, and the epoch in which the artifact occurred.

The `signal` must be a `TimeSeries` object with fields `x` (a vector of values) 
and `fs` (the sampling frequency in Hz).

Only artifacts with `test.statistic > 1000` are retained. Start and end indices 
are rounded and cast to `Int` for consistency.

Returns a DataFrame with the following columns:
- `start::Int`: Starting sample index of the artifact
- `end::Int`: Ending sample index of the artifact
- `mean.change::Float64`: Magnitude of the mean shift
- `test.statistic::Float64`: CAPA test statistic
- `epoch::Int`: Epoch number in which the artifact was found

Example:

```julia
eeg = EEG(some_eeg_file.edf)
signal = get_channel(eeg, "EEG6")
anoms = detect_artifacts(signal, 60*5)  # epoch length = 5 minutes
```
"""
function detect_artifacts(signal::TimeSeries, seg_length::Int)::Vector{Artifact}
    x = signal.x 
    fs = signal.fs
    @rput x fs seg_length
    R"result <- detect_artifacts(x, fs, seg_length)"
    @rget result
    
    start_indices = round.(Int, result[!, :start])
    end_indices = round.(Int, result[!, :end])
    mean_changes = result[!, Symbol("mean_change")]
    test_statistics = result[!, Symbol("test_statistic")]

    artifacts = Artifact[]
    for (start, stop, meanchg, teststat) in zip(start_indices, end_indices, mean_changes, test_statistics)
        epoch_s = index_to_epoch(start, fs, 30)
        epoch_e = index_to_epoch(stop, fs, 30)
        push!(artifacts, Artifact((start, stop), meanchg, teststat, (epoch_s, epoch_e)))
    end

    return artifacts
end

function detect_artifacts(eeg::EEG, channel_name::String, seg_length::Int)::Nothing
  signal = get_channel(eeg, channel_name)
  eeg._artifacts[channel_name] = detect_artifacts(signal, seg_length)
  return
end


function get_epochs_with_artifacts(artifacts::Vector{Artifact})
    epochs = Int[]
    for a in artifacts
        push!(epochs, a.epoch_range[1]:a.epoch_range[2]...)  
    end
    return sort(unique(epochs))
end


function get_epochs_with_artifacts(eeg::EEG, channel_name::String)
    artifacts = get_artifacts(eeg, channel_name)
    epochs = Int[]
    for a in artifacts
        push!(epochs, a.epoch_range[1]:a.epoch_range[2]...)  
    end
    return sort(unique(epochs))
end



"""
Incomplete but similar to previous method.
"""
function plot_artifacts_in_epochs(eeg::EEG, channel_name::String, from::Integer, to::Integer; 
                                  annotate::Bool=false)

    # Vector of artifacts in the desired region
    artifacts = filter(
      x -> x.epoch_range[1] >= from && x.epoch_range[2] <= to, 
      get_artifacts(eeg, channel_name)
    )
    if isempty(artifacts)
        @warn "No artifacts found in epoch range"
        return nothing
    end
    
    signal = get_channel(eeg, channel_name)
    s, e = epochs_to_indexes(from, to, signal.fs, 30)

    t = (s:e) ./ signal.fs
    sig = signal.x[s:e]

    # Plot base signal without legend and without artifacts
    p = plot(
        t, sig,
        label="",  # no label
        color=:black,
        xlabel="Time (s)",
        ylabel="Amplitude",
        title="Artifacts in Epochs ($from, $to)",
        legend=false,
    )

    # Overlay all artifact regions
    for art in artifacts
        s = art.loc[1]
        e = art.loc[2]

        t_anom = (s:e) ./ signal.fs
        sig_anom = signal.x[s:e]
        plot!(p, t_anom, sig_anom, lw=2, color=:red)

        # Place annotation at the middle of the artifact
        t_mid = mean(t_anom)
        y_mid = maximum(sig_anom)
        
        if annotate
          txt = "Anom \nΔμ=$(round(art.mean_change, digits=1))\nT=$(round(art.test_statistic, digits=1))"
          annotate!(p, t_mid, y_mid, text(txt, :left, 8, :red))
        end
    end

    return p
end

plot_artifacts_in_epochs(eeg, "EEG6", 200, 210)

