using FFTW
using DSP
using Statistics

"""
Structure for amplitude spectrum estimations. Estimations are by default 
one sided, with frequencies ranging from [0, fₛ/2].

The formula used is 

```math 
\\frac{2|H(f)|}{\\sum_i w_i}
```

with ``w_i`` a Hanning window.

# Fields

- `freq::Vector{<:AbstractFloat}`: Frequency range of the spectrum
- `spectrum::Vector{<:AbstractFloat}`: Estimated spectral amplitude
- `formula::String`: A string representation of the formula used for the estimation.

# Constructors
`AmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer)` : Computes a direct PSD over a signal `x` with a given `sampling_rate`.
"""
struct AmplitudeSpectrum

    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}
    formula::String

    function AmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer=0)
        N = length(x)
        hann = hanning(N) # Hanning window
        x = x .* hann
        if pad > 0
            x = zero_pad(x, pad)
        end
        ft = abs.(fft(x))
        ft = ft[1:(div(N, 2)+1)] # Make one sided
        freq = [i for i in 0:(length(ft)-1)] .* sampling_rate / N
        normalization = 2 / (sum(hann .^ 2))
        spectrum = ft * normalization
        new(freq, spectrum, "2|H(f)| / ∑ wᵢ²  with  w₁, …, wₗ a Hanning window")
    end
end


"""
Structure for PSD estimations. Estimations are by default 
one sided, with frequencies ranging from [0, fₛ/2].

The default formula is 

```math 
\\frac{2|H(f)|^2}{f_s \\sum_i w_i^2}
```

with ``w_i`` a Hanning window. This means the estimation is normalized by 
the sampling rate by default. This can be changed by setting the normalization 
parameter equal to ``\\frac{1}{f_s}``, canceling out the factor in the denominator. 

If Barlett or Welch's mehtod is used (i.e. if the second constructor is used), the formula 
becomes 


```math 
\\frac{1}{M \\varphi} \\sum_i^M \\left[ \\frac{2|H_i(f)|^2}{f_s \\sum_i w_i^2} \\right]
```

where ``w_1, \\ldots, w_n`` a Hanning window, ``M`` the number of segments,
``H_i(f)`` the FFT of the ``i``th segment of the signal, and ``\\varphi`` an
optional normalization factor defined by the `normalization` parameter
(defaults to `2 * seg_length`). Thus, with default keyword arguments, and
averaging across ``M`` windows with ``k`` samples each, the estimation is

```math 
\\frac{1}{M\times 2k} \\sum_i^M \\left[ \\frac{2|H_i(f)|^2}{f_s \\sum_i w_i^2} \\right]
```

To avoid any normalization, simply let ``\varphi = 1``.

# Fields
- `freq::Vector{<:AbstractFloat}`: Frequency range of the spectrum
- `spectrum::Vector{<:AbstractFloat}` : Estimated spectral density in dB.
- `method::String`: Estimation method used 
- `formula::String` : A string representation of the formula used for the estimation.

# Constructors
- `PSD(x::Vector, sampling_rate::Integer, pad::Integer = 0)`: Computes a direct PSD over a signal `x` with a given `sampling_rate`.
- `PSD(x::Vector, fs::Int, L::Int, overlap::Union{ <:AbstractFloat, Integer }, normalization::Union{ <:AbstractFloat, Integer } = 1)`: Splits the signal `x` into segments of length L with an `overlap` in [0, 1). The overlap is understood to be a fraction of the segment length. PSD is estimated within and averaged across all segments. If `overlap` is zero, this results in Barlett's method. If `overlap` is greater than zero, this results in Welch's method. If `pad` is zero no zero-padding is done. If `pad` is greater than zero, each segment is zero-padded to a  length of `pad`. 
- `PSD(ts::TimeSeries)`: Splits the signal `x` into segments of length L with an `overlap` in [0, 1). The overlap is understood to be a fraction of the segment length. PSD is estimated within and averaged across all segments. If `overlap` is zero, this results in Barlett's method. If `overlap` is greater than zero, this results in Welch's method. If `pad` is zero no zero-padding is done. If `pad` is greater than zero, each segment is zero-padded to a  length of `pad`. 
"""
struct PSD
    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}
    method::String
    formula::String

    function PSD(x::Vector{<:AbstractFloat}, fs::Integer; pad::Integer=0, norm_factor=1)
        N = length(x)
        hann = hanning(N) # Hanning window
        x = x .* hann
        if pad > 0
            x = zero_pad(x, pad)
        end
        ft = abs2.(fft(x))
        ft = ft[1:(div(N, 2)+1)] # Make one sided
        freq = [i for i in 0:(length(ft)-1)] .* fs / N
        normalization = 1 / (sum(hann .^ 2) * norm_factor)
        spectrum = 2 * ft * normalization
        new(freq, spectrum, "Direct (no segmentation)", "2|H(f)|² / ( ∑ wᵢ² * fₛ )  with  w₁, …, wₗ a Hanning window")
    end

    function PSD(x::Vector{<:AbstractFloat}, fs::Int, seg_length::Integer;
                overlap::Union{<:AbstractFloat,Integer} = 0.5, 
                normalization::Union{<:AbstractFloat,Integer}=-1,
                inner_normalization::Union{<:AbstractFloat, Integer}=1,
                pad::Integer=0)

        method = overlap > 0 ? "Welch's method" : "Barlett's method"
        formula = "1/(M * normalization) ∑ ᵢᴹ [ 2|Hᵢ(f)|² / ( fₛ ∑  wᵢ² ) ]  where w₁, …, wₗ a Hanning window, M the number of segments, and Hᵢ(f) the FFT of the ith segment of the signal. "

        segs = segment(x, seg_length; overlap=overlap, symmetric=true)
        M = length(segs)
        psds = map(x -> PSD(x, fs; norm_factor=inner_normalization), segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        # Default to Hans normalization: denominator = 2 * M * length(segs[1])
        if (normalization == -1)
            normalization = 2 * seg_length # seg_length = length(segs[1])
        end
        w = sum(spectrums) / (M * normalization)   
        new(freq, w, method, formula)
    end

    function PSD(ts::TimeSeries; kargs...)
        PSD(ts.x, ts.fs; kargs...)
    end
    
    function PSD(ts::TimeSeries, seg_length::Integer; kargs...)
        PSD(ts.x, ts.fs, seg_length; kargs...)
    end
end

"""
Structure for spectrogram estimation. Estimations are by default one-sided,
with frequencies ranging from [0, fₛ/2]. The signal is split into possibly overlapping 
windows of length L; within each window, a PSD method is used to compute the 
PSD with overlapping windows. 

The spectrogram is a matrix ``S^{M \\times F}`` where ``M`` is the number of windows and 
``F`` is the length of the spectrum vector in any given window (i.e. the number of 
frequencies or the frequency resolution). 

Let ``f_1, f_2, \\ldots, f_k`` be a strictly increasing sequence of
frequencies. Assume these frequencies correspond to the column indexes
``c_1, c2, \\ldots, c_k`` of ``S``. Then the mean power in the frequency range 
``[f_1, f_k]`` is

```math
\\frac{1}{M} \\sum_{i=1}^{M}\\left[\\frac{1}{c_k - c_1}\\sum_{j=c_1}^{c_k} S_{ij}\\right] = \\frac{1}{M\\big(c_k - c_1\\big)}\\sum_{i=1}^{M}\\sum_{j=c_1}^{c_k} S_{ij}
```

In this package, mean power in a frequency range is computed with the `mean_band_power` function.


# Fields
- `time::Vector` : Time domain (x)
- `freq::Vector{<:AbstractFloat}`: Frequency domain (y)
- `spectrums::Matrix{<:AbstractFloat}`: Power spectrum (z). Rows are time and columns are frequency.
- `segment_length::Integer` : Length of each segment in time.

# Constructors
- `Spectrogram(signal::Vector{<:AbstractFloat}, fs::Integer, window_length::Integer, psd_function::Function; overlap::AbstractFloat = 0.5)`: Compute the spectrogram by splitting the signal into 
potentially overlapping windows of length `window_length`. Within each window, a `psd_function` is used to compute the PSD. This function must return a `PSD` object.
`psd_overlap` parameters. An optional normalization factor and a zero-padding 
length can be included, as in the `PSD` constructor.
- `function Spectrogram(segs::Vector{Vector}, fs::Integer, segment_length::Integer, psd_function::Function; overlap::AbstractFloat = 0.5)`: This constructor does not take a signal but *a split or windowed signal*. This is useful, for example, when the EEG is split into windows corresponding to a specific period (e.g. NREM epochs). In this case, each time-instance in the Spectrogram corresponds to one of the windows. Aside from this, there is no difference with the previous constructor.
- `function Spectrogram(ts::TimeSeries, window_length::Integer, psd_function::Function; kargs...)`: Simpler constructor using a `TimeSeries` object.
"""
struct Spectrogram

    time::Vector
    freq::Vector{<:AbstractFloat}
    spectrums::Matrix{<:AbstractFloat}
    segment_length::Integer

    function Spectrogram(signal::Vector{<:AbstractFloat}, fs::Integer, window_length::Integer, psd_function::Function; 
        overlap::AbstractFloat = 0.5,
        )

        if Base.return_types(psd_function) != [PSD]
            throw(ArgumentError("The `psd_function` should be a function with return type `PSD`"))
        end

        segs = segment(signal, window_length; overlap=overlap)
        psds = map(psd_function, segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        spectrogram_data = zeros(length(spectrums), length(freq))
        for i in 1:length(spectrums)
            spectrogram_data[i, :] = spectrums[i]
        end

        new(1:length(spectrums), freq, spectrogram_data, window_length)
    end

    # If a signal is not given, but a vector of segments (e.g. a list of epochs).
    function Spectrogram(segs::Vector{Vector{T}}, psd_function::Function) where {T<:AbstractFloat}

        psds = map(psd_function, segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        spectrogram_data = zeros(length(spectrums), length(freq))
        for i in 1:length(spectrums)
            spectrogram_data[i, :] = spectrums[i]
        end

        new(1:length(spectrums), freq, spectrogram_data, length(segs[1]))
    end
    
    function Spectrogram(ts::TimeSeries, window_length::Integer, psd_function::Function; kargs...)
        Spectrogram(ts.x, ts.fs, window_length, psd_function; kargs...)
    end
end


"""
`plot_spectrogram(spec::Spectrogram, freq_lim::AbstractFloat=30.0, type::Int=1, color=:nipy_spectral)`

Plots a spectogram `spec` either in 2d (`type = 1`) or 3d (`type = 2`). An optional 
frequency limit (`freq_lim`) may be set (defaults to 30Hz). The color palette 
`color` may be set; defaults to `nipy_spectral`.
"""
function plot_spectrogram(spec::Spectrogram; freq_lim::AbstractFloat=30.0, type::Int=1, color=:nipy_spectral)

    if type == 1
        return (heatmap(spec.time, spec.freq, spec.spectrums', ylims=(0, freq_lim), color=color))
    end

    if type == 2
        return (surface(spec.time, spec.freq, spec.spectrums', ylims=(0, freq_lim),
            xlabel="Time", ylabel="Frequency (Hz)", zlabel="PSD (dB)", color=color))
    end
    throw(ArgumentError("The plot `type` argument must be either 1 (for heatmap) or 2 (for a surface plot)."))
end

"""
`plot_psd(psd::PSD; freq_lim=30.0)`

Plot a PSD with x-axis being frequency and y-axis being estimated power spectrum.
"""
function plot_psd(psd::PSD; freq_lim=30.0)
    plot(psd.freq, psd.spectrum)
end


"""
`next_power_of_two(n::Int)`

Given an integer `n`, finds the least m = 2ᵏ s.t. m ≥ n.
"""
function next_power_of_two(n::Int)
    p = 1
    while p < n
        p <<= 1
    end
    return p
end

"""
`zero_pad(v::Vector{<:AbstractFloat}, desired_length::Integer) where {T}`

Zero-pads a numeric vector `v` to a `desired_length`
"""
function zero_pad(v::Vector{T}, desired_length::Integer) where {T<:AbstractFloat}
    current_length = length(v)
    if current_length == desired_length
        return v
    end
    if desired_length < current_length
        throw(ArgumentError("Cannot zero-pad to a length inferior to the vector's length"))
    end
    padded_vector = zeros(T, desired_length)
    padded_vector[1:current_length] = v
    return padded_vector

end


"""
`freq_band(spec::Union{PSD,AmplitudeSpectrum}, lower::AbstractFloat, upper::AbstractFloat)`

Given a PSD, returns a Vector{AbstractFloat} with the powers within the frequency band [lower, upper].
"""
function freq_band(spec::Union{PSD,AmplitudeSpectrum}, lower::AbstractFloat, upper::AbstractFloat)
    indexes = findall(x -> x >= lower && x <= upper, spec.freq)
    spec.spectrum[indexes]
end


"""
`freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat, window::Integer)`

Given a spectrogram, returns a Vector{<:AbstractFloat} with the powers within a frequency band [lower, upper]
of a specific window (row of the spectrogram).
"""
function freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat, window::Integer)
    spectrum = spec.spectrums[window, :]
    indexes = findall(x -> x >= lower && x <= upper, spec.freq)
    spectrum[indexes]
end

"""
`freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)`

Given a spectrogram, returns a Matrix{<:AbstractFloat} with the powers within a frequency band [lower, upper]
across all windows.
"""
function freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)
    spectrum = spec.spectrums
    indexes = findall(x -> x >= lower && x <= upper, spec.freq)
    spectrum[:, indexes]
end

"""
`freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat, window::Integer)`

Given a spectrogram, returns the mean power in a given frequency band [lower, upper]. This function 
effectively computes 

```math
\\frac{1}{M\\big(c_k - c_1\\big)}\\sum_{i=1}^{M}\\sum_{j=c_1}^{c_k} S_{ij}
```
"""
function mean_band_power(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)
    band = freq_band(spec, lower, upper)
    # Sum columns in band and average by number of columns; 
    # i.e. compute a vector with a mean frequency power per time window
    frequency_average = mean(band, dims=2) 
    # Get the mean across time
    mean(frequency_average)
end

