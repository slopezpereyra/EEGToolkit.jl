module EEGToolkit

using Dates
using EDF
using Plots
using FFTW
using DSP
using Statistics
#using DataFrames

export Staging
export plot_hypnogram
export plot_transition_matrix

export resample 

# Artifacts (New Methods)
export buckelmueller_artifacts
export hjorth_artifacts
export hjorth_parameters

# Time Series
export TimeSeries
export segment
export epoch 
export plot_ts
export seconds_to_time
export length_in_secs
export length_in_mins
export length_in_hours
export gen_time_domain

# EEG & Masking
export EEG
export get_channel 
export get_channels 
export filter_channels 
export filter_channels!
export remove_channel! 
export plot_eeg
export get_masks
export add_mask!
export combined_mask

# PSD & Spectral Analysis
export AmplitudeSpectrum
export PSD 
export Spectrogram
export spectrum            
export plot_spectrogram 
export plot_psd
export next_power_of_two 
export zero_pad 
export freq_band 
export mean_band_power 
export total_band_power 
export relative_band_power  

# Slow Wave Detection
export SlowWave
export detect_slow_waves_massimini
export compute_morphology_metrics 
export plot_single_wave           
export plot_average_morphology   

# NREM (Assuming in nrem.jl)
export nrem

# Spindles (Assuming in spindles.jl)
export sigma_index
export relative_spindle_power

# Connectivity
export compute_wpli
export compute_coherence

# --- Includes ---

include("staging.jl")
include("hypnograms.jl")
include("ts.jl")

# Hjorth must be included before artifacts because artifacts.jl uses it
include("hjorth.jl")      
include("artifacts.jl")

include("eeg.jl")
include("resampling.jl")
include("psd.jl")
include("spindles.jl")
include("sw_detection.jl")
include("nrem.jl")
include("connectivity.jl")

end
