

"""
Structure for PSD estimations. Estimations are by default 
one sided, with frequencies ranging from [0, fₛ/2].

The default formula is 

```math 
\\frac{2|H(f)|^2}{\\zeta \\sum_i w_i^2}
```

with ``w_i`` a Hanning window and ``\\zeta`` a normalization factor which defaults 
to ``1``.  

Barlett or Welch's mehtod can be used, where the formula 
becomes 


```math 
\\frac{1}{M \\varphi} \\sum_i^M \\left[ \\frac{2|H_i(f)|^2}{ \\sum_i w_i^2} \\right]
```

where ``w_1, \\ldots, w_n`` a Hanning window, ``M`` the number of segments,
``H_i(f)`` the FFT of the ``i``th segment of the signal, and ``\\varphi``
a normalization factor defaulting to `2 * seg_length`.


# Fields
- `freq::Vector{<:AbstractFloat}`: Frequency range of the spectrum
- `spectrum::Vector{<:AbstractFloat}` : Estimated spectral density in dB.

# Constructors
- `PSD(x::Vector{<:AbstractFloat}, fs::Integer; pad::Integer=0, norm_factor=1, dB=false)`: Computes a direct PSD over a signal `x` with sampling rate `fs`. The signal may be padded to an optional length `pad`. An optional normalization factor `norm_factor` may be used. Set `dB` to true to transform the spectrum to decibels.
- `PSD(x::Vector, fs::Int, seg_length::Int; overlap::Union{ <:AbstractFloat, Integer }=0.5, normalization::Union{ <:AbstractFloat, Integer } = -1, pad::Integer=0, dB=false)`: Splits the signal `x` into segments of length `seg_length` with an `overlap` in [0, 1) (defaults to 0.5). The overlap is understood to be a fraction of the segment length. PSD is estimated within and averaged across all segments. The estimation is normalized with a `normalization` that defaults to ``2Mk'``, where `M` is the number of segments and `k'` is the number of samples in each segment (i.e. `seg_length`). Setting `overlap` to zero equates to using Barlett's method. Setting `overlap` greater than zero equates to using Welch's method. 
- `PSD(ts::TimeSeries; kargs...)` : Wrapper to apply the first constructor to a TimeSeries signal.
- `PSD(ts::TimeSeries, seg_length::Integer; kargs...)`: Wrapper to apply the second constructor (Welch or Barlett's method) to a TimeSeries signal.
"""
struct PSD
    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}

    function PSD(x::Vector{<:AbstractFloat}, fs::Integer; pad::Integer=0, norm_factor=1, dB=false)
        if pad > 0
            x = zero_pad(x, pad)
        end
        N = length(x)
        hann = hanning(N) # Hanning window
        x = x .* hann
        print(length(x))
        ft = abs2.(fft(x))
        ft = ft[1:(div(N, 2)+1)] # Make one sided
        freq = [i for i in 0:(length(ft)-1)] .* fs / N
        normalization = 1 / (sum(hann .^ 2) * norm_factor)
        spectrum = 2 * ft * normalization
        if dB 
            spectrum = pow2db.(spectrum)
        end
        new(freq, spectrum)
    end

    function PSD(x::Vector{<:AbstractFloat}, fs::Int, seg_length::Integer;
                overlap::Union{<:AbstractFloat,Integer} = 0.5, 
                normalization::Union{<:AbstractFloat,Integer}=-1,
                pad::Integer=0,
                dB = false)


        segs = segment(x, seg_length; overlap=overlap, symmetric=true)
        M = length(segs)
        psds = map(x -> PSD(x, fs), segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        # Default to Hans normalization: denominator = 2 * M * length(segs[1])
        if (normalization == -1)
            normalization = 2 * seg_length # seg_length = length(segs[1])
        end
        w = sum(spectrums) / (M * normalization)   
        if dB 
            w = pow2db.(w)
        end
        new(freq, w)
    end

    function PSD(ts::TimeSeries; kargs...)
        PSD(ts.x, ts.fs; kargs...)
    end
    
    function PSD(ts::TimeSeries, seg_length::Integer; kargs...)
        PSD(ts.x, ts.fs, seg_length; kargs...)
    end
end

"""
A spectrogram is a matrix ``S^{M \\times F}`` where ``M`` is the number of windows in 
the windowing of a signal and ``F`` is the length of the spectrum vector in any given window (i.e. the frequency resolution). It is useful to observe spectral changes in time or to compute the 
spectrum of time-regions of interest (e.g. only NREM periods in a sleep EEG). The information 
available in direct PSD can be inferred from the spectrogram with ease.

For instance, let ``f_1, f_2, \\ldots, f_k`` be a strictly increasing sequence of
frequencies. Assume these frequencies correspond to the column indexes
``c_1, c_2, \\ldots, c_k`` of ``S``. Then the mean power in the frequency range 
``[f_1, f_k]`` is

```math
\\frac{1}{M} \\sum_{i=1}^{M}\\left[\\frac{1}{c_k - c_1}\\sum_{j=c_1}^{c_k} S_{ij}\\right] = \\frac{1}{M\\big(c_k - c_1\\big)}\\sum_{i=1}^{M}\\sum_{j=c_1}^{c_k} S_{ij}
```

In this package, mean power in a frequency range is computed with the `mean_band_power` function.


# Fields
- `time::Vector` : Time domain 
- `freq::Vector{<:AbstractFloat}`: Frequency domain 
- `spectrums::Matrix{<:AbstractFloat}`: Power spectrum. Rows are time and columns are frequency; the value in `spectrums[row, freq]` is the power at time window `row` for frequency `freq`.
- `segment_length::Integer` : Length of each segment in time.

# Constructors
- `Spectrogram(segs::Vector{Vector{T}}, psd_function::Function; dB = false) where {T<:AbstractFloat}`: Given a sequence of windows ``w_1, \\ldots, w_k`` contained in the `segs` argument, computes the PSD within each window using a custom `psd_function`. 
- `Spectrogram(signal::Vector{<:AbstractFloat}, window_length::Integer, psd_function::Function; overlap::Union{AbstractFloat, Integer}=0, dB=false)`: Splits a signal into (potentially overlapping) segments of length `window_length` and computes the `Spectrogram` over this windowing using the first constructor. A custom `psd_function` is used within each window. Symmetry is enforced over the split signal, meaning that if the last segment is of length not equal to the rest, it is dropped. Thus, all windows are of equal length.
- `function Spectrogram(ts::TimeSeries, window_length::Integer, psd_function::Function; kargs...)`: Wrapper constructor for a `TimeSeries` object.
"""
struct Spectrogram
    time::Vector
    freq::Vector{<:AbstractFloat}
    spectrums::Matrix{<:AbstractFloat}
    segment_length::Integer

    function Spectrogram(segs::Vector{Vector{T}}, psd_function::Function; dB = false) where {T<:AbstractFloat}

        if Base.return_types(psd_function) != [PSD]
            throw(ArgumentError("The `psd_function` should be a function with return type `PSD`"))
        end

        psds = map(psd_function, segs)
        freq = psds[1].freq
        spectrums = [psd.spectrum for psd in psds]

        spectrogram_data = zeros(length(spectrums), length(freq))
        for i in 1:length(spectrums)
            spectrogram_data[i, :] = spectrums[i]
        end
        
        if dB 
            spectrogram_data = pow2db.(spectrogram_data)
        end

        new(1:length(spectrums), freq, spectrogram_data, length(segs[1]))
    end

    function Spectrogram(signal::Vector{<:AbstractFloat}, window_length::Integer, psd_function::Function; overlap::Union{<:AbstractFloat, Integer}=0, dB=false)
        segs = segment(signal, window_length; overlap=overlap, symmetric=true)
        Spectrogram(segs, psd_function; dB=dB)
    end
    
    function Spectrogram(ts::TimeSeries, window_length::Integer, psd_function::Function; kargs...)
        Spectrogram(ts.x, ts.fs, window_length, psd_function; kargs...)
    end
end


"""
`plot_spectrogram(spec::Spectrogram; freq_lim::AbstractFloat=30.0, type::Int=1, color=:nipy_spectral)`

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
    xlims!(0, freq_lim)
end


"""
Given an integer ``n``, finds the least ``m = 2^k`` s.t. ``m \\geq n``. 
"""
function next_power_of_two(n::Int)
    p = 1
    while p < n
        p <<= 1
    end
    return p
end

"""
`zero_pad(v::Vector{T}, desired_length::Integer) where {T<:AbstractFloat}`

Zero-pads a numeric vector `v` to a `desired_length`
"""
function zero_pad(v::Vector{T}, desired_length::Integer) where {T<:AbstractFloat}
    current_length = length(v)
    if current_length == desired_length
        return v
    end
    if desired_length < current_length
        throw(ArgumentError("Cannot zero-pad to a length inferior to the original vector's length"))
    end
    padded_vector = zeros(T, desired_length)
    padded_vector[1:current_length] = v
    return padded_vector

end


"""
`freq_band(spec::Union{PSD}, lower::AbstractFloat, upper::AbstractFloat)`

Given a `PSD`, returns a `Vector{<:AbstractFloat}` with the powers within the frequency band `[lower, upper]`.
"""
function freq_band(spec::Union{PSD}, lower::AbstractFloat, upper::AbstractFloat)
    indexes = findall(x -> x >= lower && x <= upper, spec.freq)
    spec.spectrum[indexes]
end


"""
`freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat, window::Integer)`

Given a `Spectrogram`, returns a `Vector{<:AbstractFloat}` with the powers within a frequency band `[lower, upper]`
of a specific window (row of the spectrogram).
"""
function freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat, window::Integer)
    spectrum = spec.spectrums[window, :]
    indexes = findall(x -> x >= lower && x <= upper, spec.freq)
    spectrum[indexes]
end

"""
`freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)`

Given a `Spectrogram`, returns a `Matrix{<:AbstractFloat}` with the powers within a frequency band [lower, upper]
across all time windows.
"""
function freq_band(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)
    spectrum = spec.spectrums
    indexes = findall(x -> x >= lower && x <= upper, spec.freq)
    spectrum[:, indexes]
end

"""
`mean_band_power(spec::Spectrogram, lower::AbstractFloat, upper::AbstractFloat)`

Given a `Spectrogram`, returns the mean power in a given frequency band `[lower, upper]`. This function 
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

"""
`mean_band_power(spec::PSD, lower::AbstractFloat, upper::AbstractFloat)`

Given a `PSD`, returns the mean power in a given frequency band `[lower, upper]`. 
"""
function mean_band_power(spec::PSD, lower::AbstractFloat, upper::AbstractFloat)
    band = freq_band(spec, lower, upper)
    mean(band)
end

