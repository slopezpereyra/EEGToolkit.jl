module EEGToolkit

using Dates
using EDF
using Plots
using FFTW
using DSP
using Statistics
using DataFrames
using Requires 

function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
        include("RInterface.jl")
        const EEGToolkitR = EEGToolkit.EEGToolkitR
        export EEGToolkitR
    end
end

export resample 

export Artifact 
export ArtifactData

export TimeSeries
export segment
export epoch 
export plot_ts
export seconds_to_time
export length_in_secs
export length_in_mins
export length_in_hours
export gen_time_domain

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
export total_band_power 
export analyze_eeg

export nrem

export sigma_index
export relative_spindle_power

export detect_artifacts
export plot_artifacts_in_epochs

include("ts.jl")
include("artifacts.jl")
include("eeg.jl")
include("resampling.jl")
include("psd.jl")
include("spindles.jl")
include("nrem.jl")

end
