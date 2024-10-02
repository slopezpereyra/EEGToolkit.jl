module EEGToolkit

using Dates
using EDF
using Plots
using FFTW
using DSP
using Statistics

export TimeSeries
export segment
export epoch 
export plot_ts
export seconds_to_time

export EEG
export get_channel 
export get_channels 
export filter_channels 
export filter_channels!
export remove_channel! 
export plot_eeg
export artifact_reject

export AmplitudeSpectrum
export PSD 
export Spectrogram
export plot_spectrogram 
export plot_psd
export next_power_of_two 
export zero_pad 
export freq_band 
export mean_band_power 

export nrem

export sigma_index
export relative_spindle_power

include("ts.jl")
include("eeg.jl")
include("psd.jl")
include("spindles.jl")
include("nrem.jl")

end
