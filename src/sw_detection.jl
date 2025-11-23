"""
    SlowWave

A structure representing a detected Slow Wave Oscillation.

# Fields
- `start_idx::Int`: Index of the first zero-crossing (start of negative phase).
- `neg_peak_idx::Int`: Index of the negative peak (trough).
- `mid_crossing_idx::Int`: Index of the zero-crossing between negative and positive phases.
- `pos_peak_idx::Int`: Index of the positive peak.
- `end_idx::Int`: Index of the final zero-crossing (end of positive phase).
- `neg_amp::Float64`: Amplitude of the negative peak (µV).
- `pos_amp::Float64`: Amplitude of the positive peak (µV).
- `ptp_amp::Float64`: Peak-to-peak amplitude (µV).
- `duration::Float64`: Total duration of the wave (seconds).
- `frequency::Float64`: Instantaneous frequency (Hz).
"""
struct SlowWave
    start_idx::Int
    neg_peak_idx::Int
    mid_crossing_idx::Int
    pos_peak_idx::Int
    end_idx::Int
    neg_amp::Float64
    pos_amp::Float64
    ptp_amp::Float64
    duration::Float64
    frequency::Float64
end

"""
detect_slow_waves_massimini(signal::Vector{T}, fs::Number; kwargs...)

Implements the Negative-Peak Method for Slow Wave detection as described by Massimini et al. (2004).

This algorithm anchors detection on the central negative peak (Down-state) and validates 
the wave based on zero-crossing intervals and amplitude thresholds.

# Arguments
- `signal::Vector{T}`: The EEG data vector (usually referenced to mastoids).
- `fs::Number`: Sampling frequency in Hz.

# Keyword Arguments (Thresholds)
- `freq_band::Tuple{Float64, Float64}`: Bandpass filter range. Default is `(0.1, 4.0)`.
- `amp_neg::Float64`: Min absolute amplitude for the negative trough. Default is `40.0` (µV).
                      *Note: Massimini used varying thresholds, YASA defaults to 40µV.*
- `amp_ptp::Float64`: Min peak-to-peak amplitude. Default is `75.0` (µV).
- `dur_neg::Tuple{Float64, Float64}`: Acceptable duration range for the negative half-wave. Default `(0.3, 1.5)`.
- `dur_total::Tuple{Float64, Float64}`: Acceptable duration range for the full cycle. Default `(0.5, 2.0)`.

# Returns
- `Vector{SlowWave}`: A list of detected slow wave events.

# Reference
Massimini, M., et al. (2004). "The sleep slow oscillation as a traveling wave." J Neurosci.
"""
function detect_slow_waves_massimini(
    signal::Vector{T}, 
    fs::Number;
    freq_band::Tuple{Float64, Float64}=(0.1, 4.0),
    amp_neg::Float64=40.0,
    amp_ptp::Float64=75.0,
    dur_neg::Tuple{Float64, Float64}=(0.3, 1.5),
    dur_total::Tuple{Float64, Float64}=(0.5, 2.0)
) where T <: AbstractFloat

    # 1. Filtering
    responsetype = Bandpass(freq_band[1], freq_band[2]; fs=fs)
    designmethod = Butterworth(2)
    clean_signal = filtfilt(digitalfilter(responsetype, designmethod), signal)

    # 2. Identify Zero-Crossings
    # zx_indices will contain the index *before* the crossing occurs.
    zx_indices = findall(x -> x < 0, clean_signal[1:end-1] .* clean_signal[2:end])

    detected_waves = SlowWave[]
    
    # We need at least 3 crossings to define a full cycle (Down -> Up)
    if length(zx_indices) < 3
        return detected_waves
    end

    # Iterate through crossings to find candidate negative half-waves
    i = 1
    while i <= length(zx_indices) - 2
        
        # Define candidate indices
        t1 = zx_indices[i]      # Start of potential wave
        t2 = zx_indices[i+1]    # Mid-point
        t3 = zx_indices[i+2]    # End of potential wave

        # Extract segments
        # Note: We add 1 to t2/t3 for start range to avoid double counting indices
        neg_phase_view = @view clean_signal[t1:t2]
        pos_phase_view = @view clean_signal[t2:t3]

        # CRITERION A: POLARITY CHECK
        if mean(neg_phase_view) >= 0
            i += 1
            continue
        end

        # CRITERION B: DURATION (Negative Half-Wave)
        dur_n = (t2 - t1) / fs
        if !(dur_neg[1] <= dur_n <= dur_neg[2])
            i += 1
            continue
        end

        # CRITERION C: DURATION (Full Wave)
        dur_t = (t3 - t1) / fs
        if !(dur_total[1] <= dur_t <= dur_total[2])
            i += 1
            continue
        end

        # CRITERION D: AMPLITUDE (Negative Peak)
        min_val, min_idx_rel = findmin(neg_phase_view)
        neg_peak_idx = t1 + min_idx_rel - 1

        if min_val > -abs(amp_neg) 
            i += 1
            continue
        end

        # CRITERION E: AMPLITUDE (Peak-to-Peak)
        max_val, max_idx_rel = findmax(pos_phase_view)
        pos_peak_idx = t2 + max_idx_rel - 1
        
        ptp = max_val - min_val
        
        if ptp < amp_ptp
            i += 1
            continue
        end

        # SUCCESS
        wave = SlowWave(
            t1,
            neg_peak_idx,
            t2,
            pos_peak_idx,
            t3,
            min_val,
            max_val,
            ptp,
            dur_t,
            1.0 / dur_t
        )
        push!(detected_waves, wave)

        # Move index forward. 
        i += 2
    end

    return detected_waves
end

# Wrapper for TimeSeries object
function detect_slow_waves_massimini(ts::TimeSeries; kwargs...)
    detect_slow_waves_massimini(ts.x, ts.fs; kwargs...)
end
