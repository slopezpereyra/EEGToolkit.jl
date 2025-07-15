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
    f <- function(x) { anomaly::capa(scale(x), type="mean", beta=penalty*log(length(x))) }
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


"""
    Artifact

A data structure representing a detected artifact in a time series signal.

# Fields

- `loc::Tuple{Int, Int}`: A tuple indicating the start and end sample indices of the artifact within the signal.
- `mean_change::Float64`: A numeric value representing the mean change in signal amplitude during the artifact segment. This quantifies the severity or abruptness of the artifact.
- `test_statistic::Float64`: A test statistic used to determine whether the segment should be considered an artifact.
- `epoch::Integer`: The 30-sec epoch in which the artifact begins.
- `subepoch::Integer`: The 5-sec sub-epoch (within the 30-sec epoch) where the artifact begins.

# Usage

Artifacts are typically detected using automated or manual procedures and stored in a vector of `Artifact` objects. They can be used for visualization, rejection, or correction of contaminated segments in preprocessing pipelines.
"""
struct Artifact
  loc::Tuple{Int, Int}
  mean_change::Float64 
  test_statistic::Float64
  epoch::Integer
  subepoch::Integer
end

"""
detect_artifacts(signal::TimeSeries, seg_length::Integer; penalty::Integer = 24)::Vector{Artifact}

Performs artifact detection in the `signal`. Returns a vector of `Artifact` objects.

Artifact detection is performed on each segment of the signal segmented with
`seg_length`. On each segment, the CAPA algorithm (Fisch et. al 2022) is used 
to detect epidemic distributional changes in the mean value of the segment. 
The `penalty` is an integer value β such that β ln(n) (with `n` the segment's length) 
is the penalty used by the CAPA algorithm to penalize the introduction of artifacts.
The β `penalty` defaults to `24`, which is much higher than the value recommended 
in Fisch et. al but matched human supervision on 78Hz sleep EEGs at the 
developer's laboratory. 
"""
function detect_artifacts(signal::TimeSeries, seg_length::Integer;
                          penalty::Integer = 24)::Vector{Artifact}
    x = signal.x 
    fs = signal.fs
    @rput x fs seg_length penalty
    R"result <- detect_artifacts(x, fs, seg_length)"
    @rget result
    
    start_indices = round.(Int, result[!, :start])
    end_indices = round.(Int, result[!, :end])
    mean_changes = result[!, Symbol("mean_change")]
    test_statistics = result[!, Symbol("test_statistic")]

    artifacts = Artifact[]
    for (start, stop, meanchg, teststat) in zip(start_indices, end_indices, mean_changes, test_statistics)
        epoch = index_to_epoch(start, fs, 30)
        offset = start - (epoch - 1) * (30 * fs)
        subepoch = cld(offset, 5 * fs)

        push!(artifacts, Artifact((start, stop), meanchg, teststat, epoch, subepoch))
    end

    return artifacts
end

"""
function detect_artifacts(eeg::EEG, channel_name::String, seg_length::Int)::Nothing

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
  signal = get_channel(eeg, channel_name)
  eeg._artifacts[channel_name] = detect_artifacts(signal, seg_length; penalty)
  return
end

"""
function get_epochs_with_artifacts(artifacts::Vector{Artifact})

Given a vector of artifacts, returns a vector of all 30-sec epochs which 
contain an artifact.
"""
function get_epochs_with_artifacts(artifacts::Vector{Artifact})
    epochs = Int[]
    for a in artifacts
        push!(epochs, a.epoch)  
    end
    return sort(unique(epochs))
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
function plot_artifacts_in_epochs(signal::TimeSeries, artifacts::ArtifactData, 
                                  from::Integer, to::Integer; 
                                  annotate::Bool=false)

Given a `signal` and a non-empty vector of artifacts, plots the signal and highlights
the existing artifacts from epoch `from` to epoch `to`. If `annotate` is set
to true, artifacts are annotated with their mean change and test statistic.
"""
function plot_artifacts_in_epochs(signal::TimeSeries, artifacts::ArtifactData, 
                                  from::Integer, to::Integer; 
                                  annotate::Bool=false)

    # Vector of artifacts in the desired region
    artifacts = filter(
      x -> x.epoch >= from && x.epoch <= to, 
      artifacts
    )
    if isempty(artifacts)
        @warn "No artifacts found in epoch range"
        return nothing
    end
    
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


