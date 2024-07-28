var documenterSearchIndex = {"docs":
[{"location":"#EEGToolkit.jl","page":"Home","title":"EEGToolkit.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Computational EEG analysis with emphasis in sleep neuroscience.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Developed at the Laboratory for the Study of Sleep Slow-wave activity","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Gods of the earth and sea\nSought thro' Nature to find this Tree,\nBut their search was all in vain:\nThere grows one in the Human Brain.— William Blake","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package has three aims: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Simplicity\nTransparency\nEfficiency","category":"page"},{"location":"","page":"Home","title":"Home","text":"Simplicity means that a person with little programming background should be able to use it, at least with the documentation at hand. Transparency means that any methodology implemented by the package should be accessible enough so as to be reported in a scientific paper.  Efficiency means that large EEGs (e.g. sleep EEGs) should be processed and analyzed in minutes or less.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Transparency affects primarily power spectral analysis (PSA). Most packages  don't report how PSA is done. Proprietary software is typically even more obscure. Combined with the fact that PSA is not standardized and may be computed in several many ways, this makes it very difficult to compare and rest results. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package is free software—free as in freedom. You are free to use the code as you wish and for any purpose. You are free to study the code and change it to make it do what you wish. You are free to redistribute copies of this package to help others. You are free to distribute copies of any modified version of this package. Proprietary software hinders the liberty of its users. In science, it obscures the scientific process and makes replication and collaboration difficult. If you are a scientist, use free software whenever possible.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Loading and processing EEG data\nEEG visualization\nSleep stage handling \nNREM Period detection\nPower spectral analysis\nSpindle detection algorithms","category":"page"},{"location":"#Time-series","page":"Home","title":"Time series","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TimeSeries\nsegment\nepoch\ngen_time_domain\nplot_ts\nseconds_to_time","category":"page"},{"location":"#EEGToolkit.TimeSeries","page":"Home","title":"EEGToolkit.TimeSeries","text":"A struct representing time series data.\n\nFields\n\nx::Vector{<:AbstractFloat}: Time series data.\nfs::Integer: Sampling rate.\nsecs::AbstractFloat: Duration in seconds.\nmins::AbstractFloat: Duration in minutes.\nhours::AbstractFloat: Duration in hours.\nepoch_length::Integer: Length in seconds understood to comprise an epoch (defaults to 30).\nsubepoch_length::Integer: Length in seconds understood to comprise an epoch (defaults to 30).\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.segment","page":"Home","title":"EEGToolkit.segment","text":"Splits a vector v into segments of length L with an overlap overlap expressed as a fraction of L. The overlap defaults to 0 (no overlap). Returns a vector v of vectors - i.e. Vector{Vector{T}} - with vecv_i the ith segment in the split.\n\nThe function always attempts to capture the whole vector, even if the final split is not of length L. For example, \n\n> x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]\n> segment(x, 5)\n2-element Vector{Vector{Int64}}:\n [1, 2, 3, 4, 5]\n [6, 7, 8, 9, 0]\n\n > segment(x, 7)\n2-element Vector{Vector{Int64}}:\n [1, 2, 3, 4, 5, 6, 7]\n [8, 9, 0]\n\nSet symmetric=true to ensure that, if this occurs, the last split is dropped.\n\n> segment(x, 3; symmetric=true)\n3-element Vector{Vector{Int64}}:\n [1, 2, 3]\n [4, 5, 6]\n [7, 8, 9]\n\nIf L is equal to the segment length, segment raises a warning and returns a vector with only the original vector: [v].  The return value ensures type-safety but the warning is raised because  splitting a vector over its length is potentially  a programming mistake.\n\n\n\n\n\nWrapper to segment the vector ts.x in the time  series ts.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.epoch","page":"Home","title":"EEGToolkit.epoch","text":"Returns a vector [x₁, …, xₖ] with all values xᵢ corresponding to the nth epoch in the signal.\n\n\n\n\n\nReturns a vector [x₁, …, xₖ] with all indexes corresponding to epochs n, n+1, …, m of the EEG. The default sampling rate is used to compute the indexes.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.gen_time_domain","page":"Home","title":"EEGToolkit.gen_time_domain","text":"Generates a vector of Time objects representing the time instances t_1 ldots t_n in a signal with with a given sampling rate f_s from second s to second e. For instance, f_s = 500, s = 10 e = 11 would map to 10002 10004 ldots 10998 11.\n\njulia> gen_time_domain2(500, 10, 11)\n500-element Vector{Time}:\n 00:00:10.002\n 00:00:10.004\n 00:00:10.006\n 00:00:10.008\n 00:00:10.01\n 00:00:10.012\n ⋮\n 00:00:10.992\n 00:00:10.994\n 00:00:10.996\n 00:00:10.998\n 00:00:11\n\njulia> gen_time_domain2(500, 600, 6000)\n2700000-element Vector{Time}:\n 00:10:00.002\n 00:10:00.004\n ⋮\n 01:39:59.998\n 01:40:00\n\n\n\n\n\nGenerates a vector of Time objects representing the time instances t_1 ldots t_n in a TimeSeries signal from epoch s to epoch e. \n\n\n\n\n\nGenerates a vector of Time objects representing the time instances t_1 ldots t_n in a TimeSeries signal, starting at 00:00:init where init is frac1f_s.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.plot_ts","page":"Home","title":"EEGToolkit.plot_ts","text":"Plots TimeSeries from epoch s to epoch e. The series many be normalized.\n\n\n\n\n\nPlots TimeSeries at epoch s.\n\n\n\n\n\nPlots TimeSeries. The series may be normalized.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.seconds_to_time","page":"Home","title":"EEGToolkit.seconds_to_time","text":"Helper function: maps a time in seconds to a Time object.\n\n\n\n\n\n","category":"function"},{"location":"#EEG","page":"Home","title":"EEG","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"EEG\nfilter_by_stage \nremove_channel!\nplot_eeg\nartifact_reject","category":"page"},{"location":"#EEGToolkit.EEG","page":"Home","title":"EEGToolkit.EEG","text":"A struct for the EEG data type.\n\nFields\n\nsignals::Dict{String, TimeSeries}: A dictionary mapping signal labels (strings) to arrays of floating-point values.\nstaging::Vector{String}: A vector of stage labels corresponding to each epoch.\nid::String: An identifier for the EEG.\n\nConstructors\n\nEEG(file::String; epoch_length::Integer=30, subepoch_length=5, staging::Vector{String}=[\"\"], id::String=\"\"): Instantiates an EEG from an EDF file (file).\n\nExample\n\nstaging_vector = CSV.read(\"path/to/stage_data/eeg_staging.csv\") # A vector with a stage per each epoch in the EEG\neeg_data = EEG(\"path/to/edf_data/data.edf\"; staging=staging_vector)\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.filter_by_stage","page":"Home","title":"EEGToolkit.filter_by_stage","text":"Returns all portions of an EEG channel in a given stage of the staging vector.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.remove_channel!","page":"Home","title":"EEGToolkit.remove_channel!","text":"Removes a channel from the EEG.\n\n\n\n\n\nRemoves a list of channels from the EEG.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.plot_eeg","page":"Home","title":"EEGToolkit.plot_eeg","text":"Plots EEG channels from epoch s to epoch e. Specific channels may be selected with the channels karg. The spacing argument is an added factor in the normalization of the EEG signals - the vertical distance between each signal in the plot grows proportionally to spacing.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.artifact_reject","page":"Home","title":"EEGToolkit.artifact_reject","text":"An anomaly matrix A  in  mathbbN^n  times  2  is  a  matrix  holds epoch-subepoch   pairs   that   contain    artifacts    in    a    TimeSeries. Each  row   of   (n m) of A denotes that the nth epoch contained  an artifact within sub-epoch m. \n\nThis  function  takes  a  TimeSeries  signal  and  an  anomaly  matrix  A. It   removes   from   the   signal   the   epoch-subepoch    pairs    containing artifacts. In particular, this function returns a  Vector{Vector{T}} with T<:AbstractFloat. Thus, if  result   holds   the   return   value   of    this    function,    result[i] contains   what   is   left   from    the    ith    epoch    after    removing its   contaminated   sub-epochs.   It   is   possible   that   result[i]    is empty.\n\n\n\n\n\nAn anomaly vector vecx in mathbbN^n is a sorted vector whose values are those epochs in an TimeSeries that contain  anomalies or artifacts. This function segments the TimeSeries and filters out all  epochs containing artifacts.\n\n\n\n\n\n","category":"function"},{"location":"#NREM-Period-detection","page":"Home","title":"NREM Period detection","text":"","category":"section"},{"location":"#NREM-period-definition","page":"Home","title":"NREM period definition","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Following Feinberg & Floyed and Dijk, a NREM period is a sequence of epochs satisfying the following conditions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"It starts with stages 2, 3 or 4. \nIt contains at least 15 minutes of stages 2, 3 or 4 in total.\nIt ends with 5 or more minutes of REM, or with 5 or more minutes of wakefulness. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Epochs in the sequence are allowed to contain occurrences of REM sleep or wakefulness  in between, as long as the duration of this occurrences is less than 5 minutes. But the epochs corresponding to these occurrences will not be part of the NREM period. For example, in a stage sequence of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"... - 10m of stage two - 1m of REM - 5m of stage three - 5m of REM - ...","category":"page"},{"location":"","page":"Home","title":"Home","text":"the NREM period consists of the first 10 minutes of stage 2 and the 5 minutes of stage 3, ignoring the 1 minute of REM in-between them.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Importantly, the restriction that ending REM periods must last at least 5 minutes is not imposed when detecting the first and the last NREM period in a night of sleep.","category":"page"},{"location":"#NREM-detection-algorithm","page":"Home","title":"NREM detection algorithm","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Let n be the number of epochs corresponding to 15 minutes and m the number of epochs corresponding to 5 minutes. (In 30 second epochs, n = 30 m = 10). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The algorithm assumes that the  staging field of an EEG has been set to a vector vecs that contains only the strings 1 ldots 6  (with 5 marking REM, 6 wakefulness,  unknown/unstaged).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The algorithm works by mapping vecs to alpha = s_1 ldots s_q a word over the language generated by Sigma = 1 ldots 6 .","category":"page"},{"location":"","page":"Home","title":"Home","text":"Observe that the language (5+6)^*(2+3+4)^*^* is partitioned into U and U, where U is the set of words containing at least n symbols 2 3 4 where neither 5 nor 6 occur consecutively m times. Then alpha can be decomposed into ","category":"page"},{"location":"","page":"Home","title":"Home","text":"alpha = psi_1 phi_1 psi_2 phi_2 ldots psi_k phi_k psi_k+1","category":"page"},{"location":"","page":"Home","title":"Home","text":"where phi_i = varphi_i (5^m5^* + 6^m6^*) and varphi_i in U. Such a decomposition readily provides the number of NREM periods in the EEG (i.e. k). Furthermore, the epochs which comprise these periods are easily inferable from the decomposition.","category":"page"},{"location":"","page":"Home","title":"Home","text":"nrem","category":"page"},{"location":"#EEGToolkit.nrem","page":"Home","title":"EEGToolkit.nrem","text":"Finds the k underlying NREM periods in a staging vector. Returns a vector of vectors V s.t. the ith vector in V  contains the epochs which comprise the ith NREM period. Thus, the length of V is k the number of NREM periods.\n\nThe staging field of the EEG must have been set to a vector  containing only the symbols 1, …, 6, ? where 5 denotes  REM, 6 denotes wakefulness, and ? denotes unscored/unstaged.\n\n\n\n\n\nFinds the k underlying NREM periods in the staging vector  of an EEG.\n\n\n\n\n\n","category":"function"},{"location":"#Spindle-detection","page":"Home","title":"Spindle detection","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package implements two spindle detection algorithms discussed in O'Reilly and Nielsen (2015). We give a brief overview of them here but refer to their original publications for further detail.","category":"page"},{"location":"","page":"Home","title":"Home","text":"sigma_index\nrelative_spindle_power","category":"page"},{"location":"#EEGToolkit.sigma_index","page":"Home","title":"EEGToolkit.sigma_index","text":"The sigma-index algorithm (Huupponen et al., 2007) find abnormally high amplitude values in the spindle frequency band. Per each 1 second window of the EEG, it computes\n\nthe maximum amplitude in the spindle frequency band, which we call S_max\nthe average amplitude in the low alpha and theta frequencies, which we call\n\nalpha_mean theta_mean\n\nthe maximum alpha amplitude alpha_max\n\nThe sigma-index of each window is defined to be zero if alpha_max  S_max, and otherwise\n\nf(S_max alpha_mean phi_mean) = frac2S_max alpha_mean + theta_mean  \n\nHigher values are indicative of a higher spindle probability. The rejection threshold recommended in the original paper is lambda = 45.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.relative_spindle_power","page":"Home","title":"EEGToolkit.relative_spindle_power","text":"The Relative Spindle Power (RSP) algorithm (Devuyst et al., 2011)  also detects abnormal values along the spindle frequency band.  For every 1 second window, the amplitude spectrum S(t) is computed,  and the RSP is defined as\n\nRSP(t) = fracint_11^16 S(t f) dfint_05^40 S(t f) df\n\nThis definition is more intelligible than the that of the sigma-index, insofar as it represents the ratio of the total power in the spindle band with respect to the total power in the 05 40 frequency range. It is evident by definition that 0 leq RSP leq 1. Higher values are indicative of a higher spindle probability (it should be clear that RSP is not a probability itself). The rejection threshold recommended in the original paper is lambda = 022.\n\n\n\n\n\n","category":"function"},{"location":"#Power-spectral-analysis","page":"Home","title":"Power spectral analysis","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AmplitudeSpectrum\nPSD\nplot_psd\nSpectrogram\nplot_spectrogram\nfreq_band \nmean_band_power","category":"page"},{"location":"#EEGToolkit.AmplitudeSpectrum","page":"Home","title":"EEGToolkit.AmplitudeSpectrum","text":"Structure for amplitude spectrum estimations. Estimations are by default  one sided, with frequencies ranging from [0, fₛ/2].\n\nThe formula used is \n\nfrac2H(f)sum_i w_i\n\nwith w_i a Hanning window.\n\nFields\n\nfreq::Vector{<:AbstractFloat}: Frequency range of the spectrum\nspectrum::Vector{<:AbstractFloat}: Estimated spectral amplitude\nformula::String: A string representation of the formula used for the estimation.\n\nConstructors\n\nAmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer) : Computes a direct PSD over a signal x with a given sampling_rate.\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.PSD","page":"Home","title":"EEGToolkit.PSD","text":"Structure for PSD estimations. Estimations are by default  one sided, with frequencies ranging from [0, fₛ/2].\n\nThe default formula is \n\nfrac2H(f)^2zeta sum_i w_i^2\n\nwith w_i a Hanning window and zeta a normalization factor which defaults  to 1.  \n\nBarlett or Welch's mehtod can be used, where the formula  becomes \n\nfrac1M varphi sum_i^M left frac2H_i(f)^2 sum_i w_i^2 right\n\nwhere w_1 ldots w_n a Hanning window, M the number of segments, H_i(f) the FFT of the ith segment of the signal, and varphi a normalization factor defaulting to 2 * seg_length.\n\nFields\n\nfreq::Vector{<:AbstractFloat}: Frequency range of the spectrum\nspectrum::Vector{<:AbstractFloat} : Estimated spectral density in dB.\nmethod::String: Estimation method used \nformula::String : A string representation of the formula used for the estimation.\n\nConstructors\n\nPSD(x::Vector{<:AbstractFloat}, fs::Integer; pad::Integer=0, norm_factor=1, dB=false): Computes a direct PSD over a signal x with sampling rate fs. The signal may be padded to an optional length pad. An optional normalization factor norm_factor may be used. Set dB to true to transform the spectrum to decibels.\nPSD(x::Vector, fs::Int, seg_length::Int; overlap::Union{ <:AbstractFloat, Integer }=0.5, normalization::Union{ <:AbstractFloat, Integer } = -1, pad::Integer=0, dB=false): Splits the signal x into segments of length seg_length with an overlap in [0, 1) (defaults to 0.5). The overlap is understood to be a fraction of the segment length. PSD is estimated within and averaged across all segments. The estimation is normalized with a normalization that defaults to 2Mk, where M is the number of segments and k' is the number of samples in each segment (i.e. seg_length). Setting overlap to zero equates to using Barlett's method. Setting overlap greater than zero equates to using Welch's method. \nPSD(ts::TimeSeries; kargs...) : Wrapper to apply the first constructor to a TimeSeries signal.\nPSD(ts::TimeSeries, seg_length::Integer; kargs...): Wrapper to apply the second constructor (Welch or Barlett's method) to a TimeSeries signal.\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.plot_psd","page":"Home","title":"EEGToolkit.plot_psd","text":"Plot a PSD with x-axis being frequency and y-axis being estimated power spectrum.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.Spectrogram","page":"Home","title":"EEGToolkit.Spectrogram","text":"A spectrogram is a matrix S^M times F where M is the number of windows in  the windowing of a signal and F is the length of the spectrum vector in any given window (i.e. the frequency resolution). It is useful to observe spectral changes in time or to compute the  spectrum of time-regions of interest (e.g. only NREM periods in a sleep EEG). The information  available in direct PSD can be inferred from the spectrogram with ease.\n\nFor instance, let f_1 f_2 ldots f_k be a strictly increasing sequence of frequencies. Assume these frequencies correspond to the column indexes c_1 c_2 ldots c_k of S. Then the mean power in the frequency range  f_1 f_k is\n\nfrac1M sum_i=1^Mleftfrac1c_k - c_1sum_j=c_1^c_k S_ijright = frac1Mbig(c_k - c_1big)sum_i=1^Msum_j=c_1^c_k S_ij\n\nIn this package, mean power in a frequency range is computed with the mean_band_power function.\n\nFields\n\ntime::Vector : Time domain \nfreq::Vector{<:AbstractFloat}: Frequency domain \nspectrums::Matrix{<:AbstractFloat}: Power spectrum. Rows are time and columns are frequency; the value in spectrums[row, freq] is the power at time window row for frequency freq.\nsegment_length::Integer : Length of each segment in time.\n\nConstructors\n\nSpectrogram(segs::Vector{Vector{T}}, psd_function::Function; dB = false) where {T<:AbstractFloat}: Given a sequence of windows w_1 ldots w_k contained in the segs argument, computes the PSD within each window using a custom psd_function. \nSpectrogram(signal::Vector{<:AbstractFloat}, window_length::Integer, psd_function::Function; overlap::Union{AbstractFloat, Integer}=0, dB=false): Splits a signal into (potentially overlapping) segments of length window_length and computes the Spectrogram over this windowing using the first constructor. A custom psd_function is used within each window. Symmetry is enforced over the split signal, meaning that if the last segment is of length not equal to the rest, it is dropped. Thus, all windows are of equal length.\nfunction Spectrogram(ts::TimeSeries, window_length::Integer, psd_function::Function; kargs...): Wrapper constructor for a TimeSeries object.\n\n\n\n\n\n","category":"type"},{"location":"#EEGToolkit.plot_spectrogram","page":"Home","title":"EEGToolkit.plot_spectrogram","text":"Plots a spectogram spec either in 2d (type = 1) or 3d (type = 2). An optional  frequency limit (freq_lim) may be set (defaults to 30Hz). The color palette  color may be set; defaults to nipy_spectral.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.freq_band","page":"Home","title":"EEGToolkit.freq_band","text":"Given a PSD or AmplitudeSpectrum, returns a Vector{<:AbstractFloat} with the powers within the frequency band [lower, upper].\n\n\n\n\n\nGiven a Spectrogram, returns a Vector{<:AbstractFloat} with the powers within a frequency band [lower, upper] of a specific window (row of the spectrogram).\n\n\n\n\n\nGiven a Spectrogram, returns a Matrix{<:AbstractFloat} with the powers within a frequency band [lower, upper] across all time windows.\n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.mean_band_power","page":"Home","title":"EEGToolkit.mean_band_power","text":"Given a `, returns the mean power in a given frequency band[lower, upper]`. This function  effectively computes \n\nfrac1Mbig(c_k - c_1big)sum_i=1^Msum_j=c_1^c_k S_ij\n\n\n\n\n\nGiven a PSD, returns the mean power in a given frequency band [lower, upper]. \n\n\n\n\n\n","category":"function"},{"location":"#Helpers","page":"Home","title":"Helpers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"next_power_of_two \nzero_pad ","category":"page"},{"location":"#EEGToolkit.next_power_of_two","page":"Home","title":"EEGToolkit.next_power_of_two","text":"Given an integer n, finds the least m = 2^k s.t. m geq n. \n\n\n\n\n\n","category":"function"},{"location":"#EEGToolkit.zero_pad","page":"Home","title":"EEGToolkit.zero_pad","text":"Zero-pads a numeric vector v to a desired_length\n\n\n\n\n\n","category":"function"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"#NREM-delta-power","page":"Home","title":"NREM delta power","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is an example script for computing the mean delta (delta) power in each of the NREM periods of a sleep EEG. We will use the C3 channel.","category":"page"},{"location":"","page":"Home","title":"Home","text":"# First, import the package\nusing EEGToolkit \n\n# Assuming we have the stage data in a .csv and we have some function \n# to read CSVs (e.g. from the CSV package)\nstaging = some_function_to_read_csv(\"my_staging_data.csv\")\n\n# We read an EEG that has channels C3-A2 and F3-A1. We assume the CSV had a \n# column called STAGES with the stages of each epoch.\neeg = EEG(edf_file, staging.STAGES)\n\n# We extract the TimeSeries object corresponding to C3-A2\nsignal = eeg.signals[\"C3-A2\"]\n\n# Detect the NREM periods\nnrems = nrem(eeg)\n\n# Split the C3 signal into 30-second windows (not-overlapping).\nepochs = segment(signal, signal.fs * 30)\n\n# PSD function to be used within each window in the spectrograms\npsd = x -> PSD(x, signal.fs, signal.fs * 5)\n\nmean_delta_powers = []\nfor nrem_period in nrems\n    # Extract the portion of the signal corresponding to this NREM period\n    # This is a vector of vectors [vector_1, ..., vector_k], with the ith \n    # vector being the ith epoch in this NREM period.\n    nrem_epochs = epochs[nrem_period]\n\n    # Compute spectrogram with each window being an epoch of this nrem period.\n    spec = Spectrogram(nrem_epochs, nrem_signal.fs*30, psd)\n\n    # Compute mean power in delta band (0.5 to 3.9 Hz) from the spectrogram.\n    δ = mean_band_power(spec, 0.5, 3.9)\n    # Store the result in the mean_delta_powers list.\n    push!(mean_delta_powers, δ)\nend\n\n# Now the ith element in `mean_delta_powers` is the mean delta power \n# of the ith NREM period.","category":"page"}]
}
