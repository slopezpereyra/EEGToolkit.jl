"""
    apply_bandpass(x::AbstractVector, fs::Real, freq_band::Tuple{Real, Real}; order::Integer=2)

Applies a zero-phase Butterworth band-pass filter to a vector `x`. 
Returns a new filtered vector.
"""
function apply_bandpass(x::AbstractVector, fs::Real, freq_band::Tuple{Real, Real}; order::Integer=2)
    nyquist = fs / 2
    norm_low = freq_band[1] / nyquist
    norm_high = freq_band[2] / nyquist
    
    responsetype = Bandpass(norm_low, norm_high)
    designmethod = Butterworth(order)
    return filtfilt(digitalfilter(responsetype, designmethod), x)
end

"""
    apply_bandpass!(ts::TimeSeries, freq_band::Tuple{Real, Real}; order::Integer=2)

Mutates the `ts.x` field of a `TimeSeries` object by applying a band-pass filter. 
"""
function apply_bandpass!(ts::TimeSeries, freq_band::Tuple{Real, Real}; order::Integer=2)
    ts.x .= apply_bandpass(ts.x, ts.fs, freq_band; order=order)
    return ts
end

"""
    apply_bandpass(ts::TimeSeries, freq_band::Tuple{Real, Real}; order::Integer=2)

Non-mutating version for `TimeSeries`.  Returns a new `TimeSeries` object 
containing the filtered signal.
"""
function apply_bandpass(ts::TimeSeries, freq_band::Tuple{Real, Real}; order::Integer=2)
    new_x = apply_bandpass(ts.x, ts.fs, freq_band; order=order)
    return TimeSeries(new_x, ts.fs)
end

"""
    apply_lowpass(x::AbstractVector, fs::Real, cutoff::Real; order::Integer=2)

Applies a zero-phase Butterworth low-pass filter to a vector `x`. 
Returns a new filtered vector.
"""
function apply_lowpass(x::AbstractVector, fs::Real, cutoff::Real; order::Integer=2)
    nyquist = fs / 2
    norm_cutoff = cutoff / nyquist
    responsetype = Lowpass(norm_cutoff)
    designmethod = Butterworth(order)
    return filtfilt(digitalfilter(responsetype, designmethod), x)
end

"""
    apply_lowpass!(ts::TimeSeries, cutoff::Real; order::Integer=2)

Mutates the `ts.x` field of a `TimeSeries` object by applying a low-pass filter. 
"""
function apply_lowpass!(ts::TimeSeries, cutoff::Real; order::Integer=2)
    ts.x .= apply_lowpass(ts.x, ts.fs, cutoff; order=order)
    return ts
end

"""
    apply_lowpass(ts::TimeSeries, cutoff::Real; order::Integer=2)

Non-mutating version for `TimeSeries`.  Returns a new `TimeSeries` object 
containing the filtered signal.
"""
function apply_lowpass(ts::TimeSeries, cutoff::Real; order::Integer=2)
    new_x = apply_lowpass(ts.x, ts.fs, cutoff; order=order)
    return TimeSeries(new_x, ts.fs)
end

"""
    apply_highpass(x::AbstractVector, fs::Real, cutoff::Real; order::Integer=2)

Applies a zero-phase Butterworth high-pass filter to a vector `x`. 
Returns a new filtered vector.
"""
function apply_highpass(x::AbstractVector, fs::Real, cutoff::Real; order::Integer=2)
    nyquist = fs / 2
    norm_cutoff = cutoff / nyquist
    responsetype = Highpass(norm_cutoff)
    designmethod = Butterworth(order)
    return filtfilt(digitalfilter(responsetype, designmethod), x)
end

"""
    apply_highpass!(ts::TimeSeries, cutoff::Real; order::Integer=2)

Mutates the `ts.x` field of a `TimeSeries` object by applying a high-pass filter. 
"""
function apply_highpass!(ts::TimeSeries, cutoff::Real; order::Integer=2)
    ts.x .= apply_highpass(ts.x, ts.fs, cutoff; order=order)
    return ts
end

"""
    apply_highpass(ts::TimeSeries, cutoff::Real; order::Integer=2)

Non-mutating version for `TimeSeries`.  Returns a new `TimeSeries` object.
"""
function apply_highpass(ts::TimeSeries, cutoff::Real; order::Integer=2)
    new_x = apply_highpass(ts.x, ts.fs, cutoff; order=order)
    return TimeSeries(new_x, ts.fs)
end

"""
    apply_notch(x::AbstractVector, fs::Real, freq::Real; bandwidth=2.0, order::Integer=2)

Applies a zero-phase notch (stop-band) filter to a vector `x` to attenuate a specific frequency. 
Commonly used for power-line interference (50/60 Hz).
"""
function apply_notch(x::AbstractVector, fs::Real, freq::Real; bandwidth=2.0, order::Integer=2)
    nyquist = fs / 2
    low = (freq - bandwidth/2) / nyquist
    high = (freq + bandwidth/2) / nyquist
    responsetype = Stopband(low, high)
    designmethod = Butterworth(order)
    return filtfilt(digitalfilter(responsetype, designmethod), x)
end

"""
    apply_notch!(ts::TimeSeries, freq::Real; bandwidth=2.0, order::Integer=2)

Mutates the `ts.x` field of a `TimeSeries` object by applying a notch filter at `freq`. 
"""
function apply_notch!(ts::TimeSeries, freq::Real; bandwidth=2.0, order::Integer=2)
    ts.x .= apply_notch(ts.x, ts.fs, freq; bandwidth=bandwidth, order=order)
    return ts
end

"""
    apply_notch(ts::TimeSeries, freq::Real; bandwidth=2.0, order::Integer=2)

Non-mutating version for `TimeSeries`.  Returns a new `TimeSeries` object.
"""
function apply_notch(ts::TimeSeries, freq::Real; bandwidth=2.0, order::Integer=2)
    new_x = apply_notch(ts.x, ts.fs, freq; bandwidth=bandwidth, order=order)
    return TimeSeries(new_x, ts.fs)
end
