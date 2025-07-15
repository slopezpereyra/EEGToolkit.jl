
"""
Find a rational approximation of floating point number `x` with a 
tolerance of `tol`. Returns a value of type `Rational{Integer}`.
"""
function rat_approximation(x::Float64, tol::Float64 = 1e-6)::Rational{Integer}
    max_den = floor(Int, 1 / tol)

    for denom in 1:max_den
        num = round(Int, x * denom)
        r = num // denom
        if abs(float(r) - x) < tol
            return r
        end
    end
    return round(Int, x) // 1  # fallback if nothing found
end

"""
resample(signal::TimeSeries, new_fs::Float64)::TimeSeries

Resamples the input `signal` to the given `new_fs`, using rational factor 
resampling (L/M). Returns a new `TimeSeries` with the resampled signal and 
updated sampling rate.
"""
function resample(signal::TimeSeries, new_fs::Float64)::TimeSeries
    x = signal.x
    fs = signal.fs

    # Compute rational approximation of resampling ratio
    ratio = new_fs / fs
    tol = 1e-6
    L, M = rat_approximation(ratio, tol).num, rat_approximation(ratio, tol).den

    # 1. Upsample by L (insert L-1 zeros between samples)
    upsampled = zeros(Float32, L * length(x))
    upsampled[1:L:end] .= x

    # 2. Low-pass filter at half of the lower sampling rate
    cutoff = min(fs, new_fs) / 2
    nyquist = (fs * L) / 2
    norm_cutoff = cutoff / nyquist
    filt = digitalfilter(Lowpass(norm_cutoff), FIRWindow(hamming(101)))
    filtered = filtfilt(filt, upsampled)

    # 3. Downsample by M (keep every M-th sample)
    resampled = filtered[1:M:end]

    return TimeSeries(resampled, new_fs)
end


