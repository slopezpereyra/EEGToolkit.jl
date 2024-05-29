module JEEG


export EEG
export epoch 
export time 
export plot_eeg
export plot_eeg_overlay
export get_stage_indexes
export get_stage 
export sfilter!
export artifact_reject

export overlaps 
export AmplitudeSpectrum 
export PSD 
export Spectrogram
export plot_spectrogram 
export next_power_of_two 
export zero_pad 
export freq_band 

export nrem 


include("eeg.jl")
include("nrem.jl")
include("psd.jl")

end
