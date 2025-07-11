using CSV
using DataFrames
using Plots
using RCall

# Load required R libraries and define all functions
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
df = detect_artifacts(signal, 60*5)  # epoch length = 5 minutes
```
"""
function detect_artifacts(signal::TimeSeries, epoch_length::Int)
    x = signal.x 
    fs = signal.fs
    @rput x fs epoch_length
    R"result <- detect_artifacts(x, fs, epoch_length)"
    @rget result
    return Dict(pairs(eachcol(result)))
end

function plot_artifact(i::Int, anoms::Dict, signal::TimeSeries; pad::Int=1)
    s = anoms[:start][i] 
    e = anoms[:end][i] 

    samples_per_epoch = 30 * signal.fs
    epoch_number_s = cld(s, samples_per_epoch) - pad
    epoch_number_e = cld(e, samples_per_epoch) + pad
    range = epochs_to_indexes(epoch_number_s, epoch_number_e, signal.fs, 30)
    
    plot_start = range[1]
    plot_end = range[2]

    t = (plot_start:plot_end) ./ signal.fs
    sig = signal.x[plot_start:plot_end]

    title_str = "Artifact $i (Epoch ($epoch_number_s, $epoch_number_e )"
    subtitle_str = "Δμ = $(round(anoms[:mean_change][i], digits=2)), stat = $(round(anoms[:test_statistic][i], digits=2))"
    
    p = plot(
      t, sig,
      label="Signal",
      color=:black,
      xlabel="Time (s)",
      ylabel="Amplitude",
      title=title_str * "\n" * subtitle_str
)

    # Overlay artifact region
    t_anom = (s:e) ./ signal.fs
    sig_anom = signal.x[s:e]

    plot!(t_anom, sig_anom, color=:red, label="Artifact")

    return p
end


