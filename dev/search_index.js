var documenterSearchIndex = {"docs":
[{"location":"#EEGToolkit.jl","page":"Home","title":"EEGToolkit.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Computational EEG analysis with emphasis in sleep neuroscience.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Loading EEG data\nEEG visualization\nSleep stage handling and NREM period detection\nPower spectral analysis\nSpindle detection algorithms","category":"page"},{"location":"#Function-Documentation","page":"Home","title":"Function Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"EEG\nepoch\nget_stage_indexes \nplot_eeg_overlay \nSpectrogram\ngen_time_domain\nnext_power_of_two \nget_stage \noverlaps \nAmplitudeSpectrum\nPSD\nplot_spectrogram\nzero_pad \nfreq_band \nplot_eeg \nartifact_reject","category":"page"},{"location":"#EEGToolkit.EEG","page":"Home","title":"EEGToolkit.EEG","text":"A mutable struct representing EEG data with associated metadata.\n\nFields\n\nsignals::Dict{String, Vector{<:AbstractFloat}}: A dictionary mapping signal labels (strings) to arrays of floating-point values.\nsampling_rates::Dict{String, Integer}: A dictionary mapping signals (strings) to the integer sampling rates.\nfs::Integer: A default sampling rate that will be used for calculations. Defaults to the maximum sampling rate among all signals.\nepoch_length::Integer: Length of each epoch (in seconds).\nstaging::Vector{String}: A vector of stage labels corresponding to each epoch.\nid::String: An identifier for the EEG.\n\nConstructors\n\nEEG(file, fs, epoch_length, staging): Constructs an EEG object from an EDF file (file) containing EEG data. The function reads the signals, computes the necessary metadata (fs, N, epoch_count), and initializes the EEG struct with the provided staging vector.\n\nExample\n\n```julia stagingvector = CSV.read(\"path/to/stagedata/eegstaging.csv\") # A vector with a stage per each epoch in the EEG eegdata = EEG(\"path/to/edfdata/data.edf\", 30, stagingvector)\n\nAlternatively, if no stage data exists, it is safe to do\n\neegdata = EEG(\"path/to/edfdata/data.edf\", 30, [])\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.epoch","page":"Home","title":"EEGToolkit.epoch","text":"epoch(eeg::EEG, n::Integer, fs::Integer=-1)\n\nReturns a vector [i₁, …, iₖ] with all indexes corresponding to the nth epoch of the EEG. The default sampling rate is used to compute the indexes.\n\n\n\n\n\nepoch(eeg::EEG, n::Integer, m::Integer, fs::Integer=-1)\n\nReturns a vector [i₁, …, iₖ] with all indexes corresponding to epochs n, n+1, …, m of the EEG. The default sampling rate is used to compute the indexes.\n\n\n\n\n\nepoch(eeg::EEG, n::Integer, channel::String)\n\nReturns a vector [x₁, …, xₖ] with all values of the signal channel in the nth epoch.\n\n\n\n\n\nReturns a vector [x₁, …, xₖ] with all values of the signal channel in the epochs n, n+1, …, m.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.get_stage_indexes","page":"Home","title":"EEGToolkit.get_stage_indexes","text":"get_stage_indexes(eeg::EEG, stages::Vector)\n\nThis function maps an EEG and astages vector to the array of all indexes whose values in an EEG  signal pertain to a stage in stages.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.plot_eeg_overlay","page":"Home","title":"EEGToolkit.plot_eeg_overlay","text":"plot_eeg_overlay(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)\n\nPlots with the active backend all EEG channels in channels in the  range from the nth to the mth epoch.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.Spectrogram","page":"Home","title":"EEGToolkit.Spectrogram","text":"Structure for spectrogram estimation. Estimations are by default one-sided, with frequencies ranging from [0, fₛ/2]. The signal is split into possibly overlapping  windows of length L; within each window, Welch's method is used to compute the  PSD with overlapping windows. For Barlett's method, one can set the inner window  length and the overlap to zero.\n\nFields\n\ntime::Vector : Time domain (x)\nfreq::Vector{<:AbstractFloat}: Frequency domain (y)\nspectrums::Matrix{<:AbstractFloat}: Power spectrum (z). Rows are time and columns are frequency.\nsegment_length::Integer : Length of each segment in time.\n\nConstructors\n\nSpectrogram(signal::Vector{<:AbstractFloat}, fs::Integer, segment_length::Integer, overlap::AbstractFloat, inner_window_length::Integer, inner_overlap::AbstractFloat) : Computes the Spectrogram of a signal with sampling rate fs in windows of length segment_length (in number of samples) with a certain overlap ∈ [0, 1].\n\nWithin each window, the PSD constructor is used to compute either a Welch or  a Barlett method estimation, depending on the inner_window_length and  inner_overlap parameters. An optional normalization factor and a zero-padding  length can be included, as in the PSD constructor.\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.gen_time_domain","page":"Home","title":"EEGToolkit.gen_time_domain","text":"gen_time_domain(eeg::EEG, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer}, fs::Integer=-1)\n\nGiven an EEG, generates the time vector t₁, …, tₙ corresponding to  EEG signals from time s to e.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.next_power_of_two","page":"Home","title":"EEGToolkit.next_power_of_two","text":"next_power_of_two(n::Int)\n\nGiven an integer n, finds the least m = 2ᵏ s.t. m ≥ n.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.get_stage","page":"Home","title":"EEGToolkit.get_stage","text":"function get_stage(eeg::EEG, channel::String, stages::Vector)\n\nReturns all portions of an EEG channel in a given stage of the staging vector.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.overlaps","page":"Home","title":"EEGToolkit.overlaps","text":"overlaps(v::Vector{T}, L::Int, overlap_frac::Union{Float64,Int}) where {T}\n\nSplits a vector v into segments of length L with an overlap overlap_frac expressed as a fraction of L. \n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.AmplitudeSpectrum","page":"Home","title":"EEGToolkit.AmplitudeSpectrum","text":"Structure for amplitude spectrum estimations. Estimations are by default  one sided, with frequencies ranging from [0, fₛ/2].\n\nFields\n\nfreq::Vector{<:AbstractFloat}: Frequency range of the spectrum\nspectrum::Vector{<:AbstractFloat}: Estimated spectral amplitude\nformula::String: A string representation of the formula used for the estimation.\n\nConstructors\n\nAmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer) : Computes a direct PSD over a signal x with a given sampling_rate.\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.PSD","page":"Home","title":"EEGToolkit.PSD","text":"Structure for PSD estimations. Estimations are by default  one sided, with frequencies ranging from [0, fₛ/2].\n\nThe default formula is \n\nfrac2H(f)^2f_s sum_i w_i^2\n\nwith w_i a Hanning window. This means the estimation is normalized by  the sampling rate by default. This can be changed by setting the normalization  parameter equal to rac1f_s, canceling out the factor in the denominator. \n\nIf Barlett or Welch's mehtod is used (i.e. if the second constructor is used), the formula  becomes          formula = \"1/(M * normalization) ∑ ᵢᴹ [ 2|Hᵢ(f)|² / ( fₛ ∑  wᵢ² ) ]  where w₁, …, wₗ a Hanning window, M the number of segments, and Hᵢ(f) the FFT of the ith segment of the signal. \"\n\nfrac1M varphi sum_i^M left frac2H_i(f)^2f_s sum_i w_i^2 right\n\nwhere varphi is an optional normalization factor defined by the normalization parameter (defaults to 1).\n\nFields\n\nfreq::Vector{<:AbstractFloat}: Frequency range of the spectrum\nspectrum::Vector{<:AbstractFloat} : Estimated spectral density in dB.\nmethod::String: Estimation method used \nformula::String : A string representation of the formula used for the estimation.\n\nConstructors\n\nPSD(x::Vector, sampling_rate::Integer, pad::Integer = 0): Computes a direct PSD over a signal x with a given sampling_rate.\nPSD(x::Vector, fs::Int, L::Int, overlap::Union{ <:AbstractFloat, Integer }, normalization::Union{ <:AbstractFloat, Integer } = 1): Splits the signal x into segments of length L with an overlap in [0, 1). The overlap is understood to be a fraction of the segment length. PSD is estimated within and averaged across all segments. If overlap is zero, this results in Barlett's method. If overlap is greater than zero, this results in Welch's method. If pad is zero no zero-padding is done. If pad is greater than zero, each segment is zero-padded to a  length of pad. \n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.plot_spectrogram","page":"Home","title":"EEGToolkit.plot_spectrogram","text":"plot_spectrogram(spec::Spectrogram, freq_lim::AbstractFloat=30.0, type::Int=1, color=:nipy_spectral)\n\nPlots a spectogram spec either in 2d (type = 1) or 3d (type = 2). An optional  frequency limit (freq_lim) may be set (defaults to 30Hz). The color palette  color may be set; defaults to nipy_spectral.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.zero_pad","page":"Home","title":"EEGToolkit.zero_pad","text":"zero_pad(v::Vector{<:AbstractFloat}, desired_length::Integer) where {T}\n\nZero-pads a numeric vector v to a desired_length\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.freq_band","page":"Home","title":"EEGToolkit.freq_band","text":"freq_band(spec::Union{PSD,AmplitudeSpectrum}, lower::AbstractFloat, upper::AbstractFloat)\n\nGiven a PSD, returns a Vector{AbstractFloat} with the powers within the frequency band [lower, upper].\n\n\n\n\n\nfreq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat, window::Integer)\n\nGiven a spectrogram, returns a Vector{<:AbstractFloat} with the powers within a frequency band [lower, upper] of a specific window (row of the spectrogram).\n\n\n\n\n\nfreq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)\n\nGiven a spectrogram, returns a Matrix{<:AbstractFloat} with the powers within a frequency band [lower, upper] across all windows.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.plot_eeg","page":"Home","title":"EEGToolkit.plot_eeg","text":"Plots with the active backend all EEG channels in channels in the  range from the nth to the mth epoch.\n\n\n\n\n\nplot_eeg(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)\n\nPlots with the active backend all EEG channels in channels in the  range from the nth to the mth epoch.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.artifact_reject","page":"Home","title":"EEGToolkit.artifact_reject","text":"artifact_reject(eeg::EEG, anom_matrix::Matrix, signal::String)\n\nGiven an EEG, a 2x2 matrix associating epoch-subepoch pairs with artifacts, and a signal, returns a subset of the signal with all sub-epochs containing artifacts removed.\n\nThe signal is split in epoch-length windows and each window is split in subepoch-length  windows; the matrix gives the epoch and subepoch indexes to be removed. \n\n\n\n\n\n","category":"function"}]
}