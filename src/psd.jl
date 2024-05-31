using FFTW
using DSP

"""
`overlaps(v::Vector{T}, L::Int, overlap_frac::Union{Float64,Int}) where {T}`

Splits a vector `v` into segments of length `L` with an overlap `overlap_frac` expressed as a fraction of L. 
"""
function overlaps(v::Vector{T}, L::Int, overlap_frac::Union{Float64,Int}) where {T}
    if L > length(v)
        throw(ArgumentError("Segment length L must be less than or equal to the length of the vector."))
    end

    if overlap_frac < 0 || overlap_frac >= 1.0
        throw(ArgumentError("Overlap fraction must be in the range [0, 1)."))
    end

    D = L * overlap_frac
    M = Int(ceil((length(v) - L) / (L - D)))

    segments = Vector{Vector{T}}(undef, M)
    step = Int(floor((1 - overlap_frac) * L))  # Calculate step size

    for i in 1:M
        start_idx = 1 + (i - 1) * step
        end_idx = start_idx + L - 1

        # Ensure the last segment does not exceed the length of the vector
        if end_idx > length(v)
            break
        end

        segments[i] = v[start_idx:end_idx]
    end

    return segments
end



"""
Structure for amplitude spectrum estimations. Estimations are by default 
one sided, with frequencies ranging from [0, fₛ/2].

The formula used is 

```math 
\\left[\\frac{2|H(f)|}{\\sum_i w_i}

```

with ``w_i`` a Hanning window.

# Fields

- freq::Vector{<:AbstractFloat}: Frequency range of the spectrum
- spectrum::Vector{<:AbstractFloat}: Estimated spectral amplitude
- formula::String: A string representation of the formula used for the estimation.

# Constructors
`AmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer)` : Computes a direct PSD over a signal `x` with a given `sampling_rate`.
"""
mutable struct AmplitudeSpectrum

    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}
    formula::String

    function AmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer)
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
        formula = "1/(M * normalization) ∑ ᵢᴹ [ 2|Hᵢ(f)|² / ( fₛ ∑  wᵢ² ) ]  where w₁, …, wₗ a Hanning window, M the number of segments, and Hᵢ(f) the FFT of the ith segment of the signal. "

```math 
\\frac{1}{M \\varphi} \\sum_i^M \\left[ \\frac{2|H_i(f)|^2}{f_s \\sum_i w_i^2} \\right]
```

where ``\\varphi`` is an optional normalization factor defined by the `normalization` parameter (defaults to 1).

# Fields
- `freq::Vector{<:AbstractFloat}`: Frequency range of the spectrum
- `spectrum::Vector{<:AbstractFloat}` : Estimated spectral density in dB.
- `method::String`: Estimation method used 
- `formula::String` : A string representation of the formula used for the estimation.

# Constructors
- `PSD(x::Vector, sampling_rate::Integer, pad::Integer = 0)`: Computes a direct PSD over a signal `x` with a given `sampling_rate`.
- `PSD(x::Vector, fs::Int, L::Int, overlap::Union{ <:AbstractFloat, Integer }, normalization::Union{ <:AbstractFloat, Integer } = 1)`: Splits the signal `x` into segments of length L with an `overlap` in [0, 1). The overlap is understood to be a fraction of the segment length. PSD is estimated within and averaged across all segments. If `overlap` is zero, this results in Barlett's method. If `overlap` is greater than zero, this results in Welch's method. If `pad` is zero no zero-padding is done. If `pad` is greater than zero, each segment is zero-padded to a  length of `pad`. 
"""
struct PSD
    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}
    method::String
    formula::String

    function PSD(x::Vector{<:AbstractFloat}, fs::Integer, pad::Integer=0)
        N = length(x)
        hann = hanning(N) # Hanning window
        x = x .* hann
        if pad > 0
            x = zero_pad(x, pad)
        end
        ft = abs2.(fft(x))
        ft = ft[1:(div(N, 2)+1)] # Make one sided
        freq = [i for i in 0:(length(ft)-1)] .* fs / N
        normalization = 1 / (sum(hann .^ 2) * fs)
        spectrum = 2 * ft * normalization
        spectrum = pow2db.(spectrum)
        new(freq, spectrum, "Direct (no segmentation)", "2|H(f)|² / ( ∑ wᵢ² * fₛ )  with  w₁, …, wₗ a Hanning window")
    end

    function PSD(x::Vector{<:AbstractFloat}, fs::Int, L::Int, overlap::Union{<:AbstractFloat,Integer}, normalization::Union{<:AbstractFloat,Integer}=1,
        pad::Integer=0)


        if (L >= length(x))
            throw(ArgumentError("The length `L` of each segment should be inferior 
            to the length of the vector."))
        end
        if (overlap < 0 || overlap >= 1)
            throw(ArgumentError("The `overlap` should range in (0, 1)"))
        end

        method = overlap > 0 ? "Welch's method" : "Barlett's method"
        formula = "1/(M * normalization) ∑ ᵢᴹ [ 2|Hᵢ(f)|² / ( fₛ ∑  wᵢ² ) ]  where w₁, …, wₗ a Hanning window, M the number of segments, and Hᵢ(f) the FFT of the ith segment of the signal. "

        function H(signal::Vector, window::Vector, seg_length::Integer)
            signal = signal .* window
            if pad > 0
                signal = zero_pad(signal, pad)
            end
            ft = abs2.(fft(signal))
            ft = ft[1:(div(seg_length, 2)+1)] # One sided
            return 2 .* ft ./ (sum(window .^ 2) * fs)
        end

        hann = hanning(L)
        segs = overlaps(x, L, overlap)
        M = length(segs)
        w = sum(map(x -> H(x, hann, L), segs)) ./ (M * normalization)   # Hans uses denominator 2 * M * length(segs[1])
        w = pow2db.(w)
        freq = [i for i in 0:(div(L, 2))] .* fs / L
        new(freq, w, method, formula)
    end
end

"""
Structure for spectrogram estimation. Estimations are by default one-sided,
with frequencies ranging from [0, fₛ/2]. The signal is split into possibly overlapping 
windows of length L; within each window, Welch's method is used to compute the 
PSD with overlapping windows. For Barlett's method, one can set the inner window 
length and the overlap to zero.

# Fields
- `time::Vector` : Time domain (x)
- `freq::Vector{<:AbstractFloat}`: Frequency domain (y)
- `spectrums::Matrix{<:AbstractFloat}`: Power spectrum (z). Rows are time and columns are frequency.
- `segment_length::Integer` : Length of each segment in time.

# Constructors
- `Spectrogram(signal::Vector{<:AbstractFloat}, fs::Integer, segment_length::Integer, overlap::AbstractFloat, inner_window_length::Integer, inner_overlap::AbstractFloat)` : Computes the `Spectrogram` of a `signal` with sampling rate `fs` in windows of length `segment_length` (in number of samples) with a certain `overlap` ∈ [0, 1].
Within each window, the `PSD` constructor is used to compute either a Welch or 
a Barlett method estimation, depending on the `inner_window_length` and 
`inner_overlap` parameters. An optional normalization factor and a zero-padding 
length can be included, as in the `PSD` constructor.
"""
struct Spectrogram

    time::Vector
    freq::Vector{<:AbstractFloat}
    spectrums::Matrix{<:AbstractFloat}
    segment_length::Integer

    function Spectrogram(signal::Vector{<:AbstractFloat}, fs::Integer, segment_length::Integer, overlap::AbstractFloat,
        inner_window_length::Integer, inner_overlap::AbstractFloat, normalization::Union{AbstractFloat,Integer}=1,
        pad::Integer=0)

        segs = overlaps(signal, segment_length, overlap)
        psds = map(x -> PSD(x, fs, inner_window_length, inner_overlap, normalization, pad), segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        spectrogram_data = zeros(length(spectrums), length(freq))
        for i in 1:length(spectrums)
            spectrogram_data[i, :] = spectrums[i]
        end

        new(1:length(spectrums), freq, spectrogram_data, segment_length)
    end

    # If a signal is not given, but a vector of segments (e.g. a list of epochs).
    function Spectrogram(segs, fs::Integer, segment_length::Integer, overlap::AbstractFloat,
        normalization::Union{AbstractFloat,Integer}=1,
        pad::Integer=0)

        psds = map(x -> PSD(x, fs, segment_length, overlap, normalization, pad), segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        spectrogram_data = zeros(length(spectrums), length(freq))
        for i in 1:length(spectrums)
            spectrogram_data[i, :] = spectrums[i]
        end

        new(1:length(spectrums), freq, spectrogram_data, segment_length)
    end
end


"""
`plot_spectrogram(spec::Spectrogram, freq_lim::AbstractFloat=30.0, type::Int=1, color=:nipy_spectral)`

Plots a spectogram `spec` either in 2d (`type = 1`) or 3d (`type = 2`). An optional 
frequency limit (`freq_lim`) may be set (defaults to 30Hz). The color palette 
`color` may be set; defaults to `nipy_spectral`.
"""
function plot_spectrogram(spec::Spectrogram, freq_lim::AbstractFloat=30.0, type::Int=1, color=:nipy_spectral)

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
function zero_pad(v::Vector{<:AbstractFloat}, desired_length::Integer) where {T}
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

