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
- `epoch::Integer`: To which epoch does the (beginning of the) slow wave belong? Useful for masking.
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
    epoch::Integer
end

"""
detect_slow_waves(signal::Vector{T}, fs::Number; kwargs...)

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
function detect_slow_waves(
    signal::Vector{T}, 
    fs::Number;
    freq_band::Tuple{Float64, Float64}=(0.1, 4.0),
    amp_neg::Float64=40.0,
    amp_ptp::Float64=75.0,
    dur_neg::Tuple{Float64, Float64}=(0.3, 1.5),
    dur_total::Tuple{Float64, Float64}=(0.5, 2.0)
) where T <: AbstractFloat

    clean_signal = apply_bandpass(signal, fs, freq_band)

    # 2. Identify Zero-Crossings
    zx_indices = findall(x -> x < 0, clean_signal[1:end-1] .* clean_signal[2:end])

    detected_waves = SlowWave[]
    
    if length(zx_indices) < 3
        return detected_waves
    end

    i = 1
    while i <= length(zx_indices) - 2
        t1 = zx_indices[i]
        t2 = zx_indices[i+1]
        t3 = zx_indices[i+2]

        neg_phase_view = @view clean_signal[t1:t2]
        pos_phase_view = @view clean_signal[t2:t3]

        # Criteria Checks
        if mean(neg_phase_view) >= 0; i += 1; continue; end
        
        dur_n = (t2 - t1) / fs
        if !(dur_neg[1] <= dur_n <= dur_neg[2]); i += 1; continue; end

        dur_t = (t3 - t1) / fs
        if !(dur_total[1] <= dur_t <= dur_total[2]); i += 1; continue; end

        min_val, min_idx_rel = findmin(neg_phase_view)
        if min_val > -abs(amp_neg); i += 1; continue; end

        max_val, max_idx_rel = findmax(pos_phase_view)
        ptp = max_val - min_val
        if ptp < amp_ptp; i += 1; continue; end

        # Calculate Epoch (cld = ceiling division)
        samples_per_epoch = fs * 30
        current_epoch = Int(cld(t1, samples_per_epoch))

        # Register Wave
        wave = SlowWave(
            t1, 
            t1 + min_idx_rel - 1, 
            t2, 
            t2 + max_idx_rel - 1, 
            t3,
            min_val, 
            max_val, 
            ptp, 
            dur_t, 
            1.0 / dur_t,
            current_epoch
        )
        push!(detected_waves, wave)
        i += 2
    end
    return detected_waves
end

"""
    detect_slow_waves(ts::TimeSeries; kwargs...)
Wrapper for `TimeSeries` input. Extracts the signal and sampling frequency, then calls the main detection function.
"""
function detect_slow_waves(ts::TimeSeries; kwargs...)
    detect_slow_waves(ts.x, ts.fs; kwargs...)
end



"""
    plot_single_wave(signal, fs, wave; padding=2.0, freq_band=(0.1, 4.0))

Plots a specific detected Slow Wave with context. The plot includes:
- The filtered signal segment around the wave (in gray).
- The detected wave itself highlighted in red.
- Markers for the negative and positive peaks.
- Vertical dashed lines for the zero-crossings (start, mid, end).
"""
function plot_single_wave(signal::Vector, fs::Number, wave; 
                          padding::Float64=2.0, 
                          freq_band::Tuple{Float64, Float64}=(0.1, 4.0))
    
    # 1. Calculate context window with extra buffer for filtering
    # A buffer prevents edge artifacts from the filter appearing in the plot
    filter_buffer = round(Int, 5.0 * fs) 
    
    pad_samples = round(Int, padding * fs)
    
    # Indices for the visual window
    plot_start = max(1, wave.start_idx - pad_samples)
    plot_end = min(length(signal), wave.end_idx + pad_samples)
    
    # Indices for the filtering window (Visual + Buffer)
    filt_start = max(1, plot_start - filter_buffer)
    filt_end = min(length(signal), plot_end + filter_buffer)
    
    # Extract raw segment
    raw_segment = signal[filt_start:filt_end]
    
    # 2. Apply Filter (using the helper with fixed Bandpass)
    clean_segment = apply_bandpass(raw_segment, fs, freq_band)
    
    # 3. Slice the clean segment back to the plotting window
    # Map global indices to local filtered segment indices
    local_plot_start = plot_start - filt_start + 1
    local_plot_end = plot_end - filt_start + 1
    
    sig_context = clean_segment[local_plot_start:local_plot_end]
    
    # Indices of the wave relative to the context vector
    wave_rel_start = wave.start_idx - plot_start + 1
    wave_rel_end = wave.end_idx - plot_start + 1
    
    sig_wave = sig_context[wave_rel_start:wave_rel_end]
    
    # 4. Create Time Vectors (in seconds)
    # Align time axis so the wave start is meaningful or just absolute time
    t_context = (plot_start:plot_end) ./ fs
    t_wave = (wave.start_idx:wave.end_idx) ./ fs
    
    # 5. Create Plot
    p = plot(t_context, sig_context, 
        label="Context (Filtered)", 
        color=:gray, 
        alpha=0.6, 
        lw=1.5,
        title="Slow Wave Event (Dur: $(round(wave.duration, digits=2))s)",
        xlabel="Time (s)", 
        ylabel="Amplitude (µV)",
        legend=:topright
    )
    
    # Overlay the detected wave in Red
    plot!(p, t_wave, sig_wave, 
        label="Detected SW", 
        color=:red, 
        lw=2.5
    )
    
    # Mark Peaks
    scatter!(p, [wave.neg_peak_idx/fs], [wave.neg_amp], color=:blue, ms=6, label="Neg Peak")
    scatter!(p, [wave.pos_peak_idx/fs], [wave.pos_amp], color=:green, ms=6, label="Pos Peak")
    
    # Mark Zero Crossings
    vline!(p, [wave.start_idx/fs, wave.mid_crossing_idx/fs, wave.end_idx/fs], 
        color=:black, 
        linestyle=:dash, 
        label=""
    )
    
    # Draw 0uV line
    hline!(p, [0], color=:black, lw=0.5, label="")
    
    return p
end

"""
    plot_single_wave(ts::TimeSeries, wave; kwargs...)
Wrapper for `TimeSeries` input. Extracts the signal and sampling frequency, then calls the main plotting function.
"""
function plot_single_wave(ts::TimeSeries, wave; kwargs...)
    plot_single_wave(ts.x, ts.fs, wave; kwargs...)
end


"""
    plot_average_morphology(signal::Vector, fs::Number, waves::Vector; window=1.0, freq_band=(0.1, 4.0), align=:neg_peak)

Calculates and plots the Grand Average (mean waveform) of a set of detected slow
waves, aligned to a specific feature.

This function allows you to visualize the "prototypical" shape of the slow
oscillations in your data by averaging all detected events. It handles signal
filtering, segment extraction, alignment, and visualization with standard
deviation confidence intervals.

# Arguments

- `signal::Vector`: The raw or pre-processed EEG data vector (single channel).
- `fs::Number`: The sampling frequency of the signal in Hz.
- `waves::Vector`: A list of detected wave objects (e.g., `SlowWave` structs). Each object must contain index fields corresponding to the chosen alignment (e.g., `neg_peak_idx`).

# Keyword Arguments

- `window::Float64`: The time window (in seconds) to include *before* and *after* the alignment point. 
    - Default: `1.0` (resulting in a 2-second total plot duration: -1.0s to +1.0s).
    - **Note:** Waves that are too close to the start or end of the signal to fit this window will be excluded.

- `freq_band::Tuple{Float64, Float64}`: The frequency range for the bandpass filter applied prior to averaging.
    - Default: `(0.1, 4.0)` Hz (Delta band).
    - **Purpose:** Ensures the morphology reflects the slow-wave component specifically, removing noise or faster activity (like spindles) that might obscure the mean shape.

- `align::Symbol`: The feature of the wave to use as the center (0s) of the time axis.
    - `:neg_peak` (Default): Aligns to the trough of the slow wave (the "Down-state").
    - `:pos_peak`: Aligns to the subsequent positive peak (the "Up-state").
    - `:mid`: Aligns to the zero-crossing between the negative and positive phases.
    - `:start`: Aligns to the first zero-crossing (start of the wave).

# Returns

- `Plots.Plot`: A plot object displaying the mean waveform (thick blue line) and the standard deviation (shaded region).
"""
function plot_average_morphology(signal::Vector, fs::Number, waves::Vector; 
                                 window::Float64=1.0, 
                                 freq_band::Tuple{Float64, Float64}=(0.1, 4.0),
                                 align::Symbol=:neg_peak)
    
    # 1. Filter the entire signal (or large chunks) first for efficiency
    # Note: If signal is massive, you might want to filter per segment, 
    # but for typical sleep EEG, filtering the whole vector is fine.
    println("Filtering signal for Grand Average...")
    clean_signal = apply_bandpass(signal, fs, freq_band)
    
    win_samples = round(Int, window * fs)
    # Time axis centered at 0
    time_axis = range(-window, window, length=2*win_samples+1)
    
    # Collect valid segments
    segments = Vector{Float64}[]
    
    for w in waves
        # Determine center index based on alignment choice
        center = if align == :mid
            w.mid_crossing_idx
        elseif align == :start
            w.start_idx
        elseif align == :pos_peak
            w.pos_peak_idx
        else # Default to :neg_peak
            w.neg_peak_idx
        end
        
        # Extract if bounds permit
        start_idx = center - win_samples
        end_idx = center + win_samples
        
        if start_idx >= 1 && end_idx <= length(clean_signal)
            push!(segments, clean_signal[start_idx:end_idx])
        end
    end
    
    if isempty(segments)
        @warn "No waves fit within the window bounds."
        return plot()
    end
    
    # Stack into Matrix (Rows=Time, Cols=Waves)
    # hcat expects vectors, so we get (Time x N)
    seg_matrix = hcat(segments...)
    
    # Compute Statistics across columns (dims=2)
    mean_wave = vec(mean(seg_matrix, dims=2))
    std_wave = vec(std(seg_matrix, dims=2))
    
    # Determine Labels
    title_text = "Grand Average (N=$(length(segments)), Aligned: $(align))"
    xlabel_text = "Time relative to alignment (s)"
    
    # Plotting
    p = plot(time_axis, mean_wave, 
        ribbon=std_wave, 
        fillalpha=0.2, 
        color=:blue, 
        lw=3, 
        label="Mean ± SD",
        title=title_text,
        xlabel=xlabel_text,
        ylabel="Amplitude (µV)"
    )
    
    # Add guidelines
    vline!(p, [0], color=:black, linestyle=:dash, label="Center")
    hline!(p, [0], color=:black, lw=1, label="")
    
    # Annotations
    if align == :neg_peak
        annotate!(p, 0, minimum(mean_wave), text("↓ Down-state", :top, 8))
    end
    
    return p
end

"""
    plot_average_morphology(ts::TimeSeries, waves::Vector; kwargs...)
Wrapper for `TimeSeries` input. Extracts the signal and sampling frequency, then calls the main plotting function.
"""
function plot_average_morphology(ts::TimeSeries, waves::Vector; kwargs...)
    plot_average_morphology(ts.x, ts.fs, waves; kwargs...)
end

"""
    compute_morphology_metrics(waves::Vector{SlowWave}, signal::AbstractVector, fs::Number)

Computes aggregate morphology statistics for a set of waves. 
Total duration is inferred directly from the length of the provided `signal`.

Calculates:
1. Density (Count / min)
2. Mean Duration
3. Mean PTP Amplitude
4. Mean Slope (Down-state slope)
"""
function compute_morphology_metrics(waves::Vector{SlowWave}, signal::AbstractVector, fs::Number)
    if isempty(waves)
        return (density=0.0, mean_dur=NaN, mean_ptp=NaN, mean_slope=NaN)
    end

    # DEDUCTION: Calculate duration from the signal itself
    total_duration_sec = length(signal) / fs
    duration_min = total_duration_sec / 60.0
    
    # 1. Density
    count = length(waves)
    density = count / duration_min

    # 2. Mean Duration
    mean_dur = mean([w.duration for w in waves])

    # 3. Mean PTP
    mean_ptp = mean([w.ptp_amp for w in waves])

    # 4. Mean Slope (1st Segment: Start -> NegPeak)
    slopes = Float64[]
    for w in waves
        # Calculate time delta for the descending slope
        dt = (w.neg_peak_idx - w.start_idx) / fs
        
        # Avoid division by zero (though unlikely in valid waves)
        if dt > 0
            push!(slopes, abs(w.neg_amp) / dt)
        end
    end
    mean_slope = isempty(slopes) ? NaN : mean(slopes)

    return (density, mean_dur, mean_ptp, mean_slope)
end

"""
    compute_morphology_metrics(waves::Vector{SlowWave}, ts::TimeSeries)

Wrapper for `TimeSeries` input. Extracts the signal and sampling frequency, then calls the main computation function.
"""
function compute_morphology_metrics(waves::Vector{SlowWave}, ts::TimeSeries)
    compute_morphology_metrics(waves, ts.x, ts.fs)
end

"""
    filter_waves(waves::Vector{SlowWave}, mask::BitVector)

Filters a list of SlowWaves, returning only those that belong to epochs 
where `mask` is true.
"""
function filter_waves(waves::Vector{SlowWave}, mask::Union{BitVector, Vector{Bool}})
    filter(w -> 1 <= w.epoch <= length(mask) && mask[w.epoch], waves)
end

