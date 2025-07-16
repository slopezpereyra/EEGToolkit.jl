# src/RInterface.jl
module EEGToolkitR
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
`detect_artifacts(signal::TimeSeries, seg_length::Integer; penalty::Integer = 24)::Vector{Artifact}`

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
function anomaly(x::Vector{<:AbstractFloat}, fs::Integer, seg_length::Integer;
                          penalty::Integer = 24)
    @rput x fs seg_length penalty
    R"result <- detect_artifacts(x, fs, seg_length)"
    @rget result
    
    start_indices = round.(Int, result[!, :start])
    end_indices = round.(Int, result[!, :end])
    mean_changes = result[!, Symbol("mean_change")]
    test_statistics = result[!, Symbol("test_statistic")]

    # Horrible, but necessary to keep isolation between modules. 
    # These four vectors will be used in `detect_artifacts` (eeg.jl)
    # to create artifact objects, but since these operations depend on 
    # the EEGToolkit module they cannot be performed here.
    return start_indices, end_indices, mean_changes, test_statistics
end


end
