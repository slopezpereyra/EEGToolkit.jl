# EEGToolkit

[![Build Status](https://github.com/slopezpereyra/JEEG.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/slopezpereyra/JEEG.jl/actions/workflows/CI.yml?query=branch%3Amain)


> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of
> Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

A scientific package for computational EEG analysis with an
emphasis on methodological transparency. Current features:

- EEG Computational Toolkit
    - Loading EEG data
    - EEG visualization
    - Sleep stage handling and NREM period detection
    - Power spectral analysis
    - Spindle detection algorithms

Read the documentation [here](https://slopezpereyra.github.io/EEGToolkit.jl/dev/).

### What is an EEG 

In this package, an EEG is a struct with fields 

- `signals`: A dictionary mapping channel names (`String`) to numeric signals `Vector{<:AbstractFloat}`
- `sampling_rates`: A dictionary mapping channel names (`String`) to their sampling rates (`Integer`)
- `epoch_length`: Number of seconds (`Integer`) understood to comprise an epoch.
- `staging`: A `Vector{String}` s.t. the $i$th word in the vector denotes the stage of the $i$th epoch in the EEG. 
- `id`: A `String` identifier for this EEG; defaults to an empty string.

### Reading EEG data

To read an EDF file, 

```julia
eeg_var_name = EEG(file, epoch_length, staging, id)
```

The `EEG` constructor instantiates an `EEG` structure.

### Spectral analysis

It is straightforward to compute the spectrogram of a signal.


```julia
S = Spectrogram(signal, eeg.fs, 3, 0.5) # Compute spectrogram with 3 second segments and 0.5 segment overlap.
p = plot_spectrogram(S, 30.0, 2) # Plot the spectrogram with limit frequency 30.0; type 2 plot = surface plot.
```

![Image](imgs/spetrogram_plot.png)

We can also take a look at the spectrogram as a heatmap:

```julia
p = plot_spectrogram(S, 30.0, 1, :inferno) # Color scheme inferno is better for heatmap
```

![Image](imgs/spetrogram_hplot.png) 

The power spectrum is easily computed and easily plotted. It is easy to set Welch's method, Barlett's method, 
or direct (no segmentation) PSD estimation. In this case we use Welch's method with 3 second windows and $0.5$ overlap.

```julia
psd = PSD(signal, eeg.fs, eeg.fs * 3, 0.5)
plot(psd.freq, pow2db.(psd.spectrum), xlab="Frequency (Hz)", ylab="PSD (dB)", legend=false)
```

The `psd` function has an optional last argument (left unwritten in our example) called 
`normalization`. It represents any additional normalization factor to be used in normalizing 
the PSD.

For transparency, the `PSD` struct contains the fields `method` and `formula`, which 
hold methodological information about how the PSD was computed. This is intended to avoid 
obscurities typical of other packages (for example, what normalization factors 
are being used). In this case,

```
julia> psd.method
"Welch's method"

julia> psd.formula
"1/(M * normalization) ∑ ᵢᴹ [ 2|Hᵢ(f)|² / ∑  wᵢ² ]  where w₁, …, wₗ a Hanning window, M the number of segments, and Hᵢ(f) the FFT of the ith segment of the signal. "
```


![Image](imgs/psd.png) 

### Signal filtering 

Our filtering functions are simply wrappers around the `DSP.jl` package 
and are very easy to use. The dispatches of the `filter!` function are:

```julia
sfilter!(eeg::EEG, channel::String, digfilter, cut_off)

sfilter!(eeg::EEG, channels::Vector{<:String}, digfilter, cut_off)

sfilter!(eeg::EEG, digfilter, cut_off)
```

If no channel is given, all EEG signals are filtered. For example,
`filter!(eeg, Lowpass, 1)` applies a low-pass filter with cut-off frequency
$1$Hz to all EEG signals.

### Spindle detection 








