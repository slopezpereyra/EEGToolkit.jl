"""
    epoch_window_bounds(n_epochs, epoch_index, window_length) -> (start_idx, end_idx)

Compute the start and end indices of a centered sliding window over a sequence
of epochs.

The window is centered on `epoch_index` and has total length
`window_length`. If the ideal centered window would extend beyond the
available epoch range `[1, n_epochs]`, the window is shifted to remain
within valid bounds while preserving the requested window length whenever
possible.

This function is marked with `@inline` because it is intended for use
inside tight loops where minimizing function call overhead improves
performance.

# Arguments

- `n_epochs::Int`: Total number of available epochs.
- `epoch_index::Int`: Index of the center epoch.
- `window_length::Int`: Desired window size (must be odd).

# Returns

- `(start_idx, end_idx)`: Tuple of integer indices defining the window bounds.

# Notes

- The window length must be odd to ensure a symmetric window around
  the center epoch under ideal conditions.
- Near the boundaries, the window may be shifted so that indices remain valid.

# Throws

- `ArgumentError` if `window_length` is even.
"""
@inline function epoch_window_bounds(
    n_epochs::Int,
    epoch_index::Int,
    window_length::Int
)
    isodd(window_length) || throw(ArgumentError("window_length must be odd"))

    radius = window_length ÷ 2
    start_idx = epoch_index - radius

    if start_idx < 1
        start_idx = 1
    end

    if (start_idx + window_length - 1) > n_epochs
        start_idx = max(1, n_epochs - window_length + 1)
    end

    end_idx = min(n_epochs, start_idx + window_length - 1)

    return start_idx, end_idx
end


"""
    buckelmueller_artifacts(signal::TimeSeries;
                            window_length=15,
                            delta_threshold=2.5,
                            beta_threshold=2.0) -> Vector{Bool}

Detect artifact epochs using a local spectral power ratio criterion inspired by
Buckelmueller-style artifact detection.

The signal is segmented into 30-second epochs. 
For each epoch, power spectral density (PSD) is computed and band power is
extracted for:

- Delta band: 0.6–4.6 Hz
- Beta band: 40–60 Hz

An epoch is flagged as an artifact if its band power exceeds the mean band power
of neighboring epochs within a sliding window (excluding the center epoch) by a
specified ratio threshold.

Specifically, an epoch `i` is marked as an artifact if:

    delta_power[i] / local_delta_mean > delta_threshold

or

    beta_power[i] / local_beta_mean > beta_threshold

where local means are computed over the surrounding window.


Segmentation is done symmetrically, so if the signal has an ending epoch of <30
seconds, it will will be dropped. This may result in the output mask having a
length of `num_epochs(signal) - 1` instead of `num_epochs(signal)`.


# Arguments

- `signal::TimeSeries`: Input time series containing signal samples and sampling frequency.

# Keyword Arguments

- `window_length::Int=15`: Number of epochs used to compute the local average.
  Must be odd.
- `delta_threshold::Real=2.5`: Threshold for delta-band power ratio.
- `beta_threshold::Real=2.0`: Threshold for beta-band power ratio.

# Returns

- `Vector{Bool}`: Boolean mask where `true` indicates an epoch classified as artifact.

# Implementation Notes

- PSD and band power are computed once per epoch for efficiency.
- Local averages exclude the center epoch to avoid bias.
- Uses `@inbounds` to avoid bounds checking inside inner loops.

# Assumptions

- Epoch duration is fixed at 30 seconds.
- Sampling frequency is provided by `signal.fs`.
"""
function buckelmueller_artifacts(
    signal::TimeSeries;
    window_length::Int = 15,
    delta_threshold::Real = 2.5,
    beta_threshold::Real = 2.0
)

    # Segment signal into 30s epochs
    fs = signal.fs
    epoch_len = fs * 30

    segments = segment(signal.x, epoch_len; symmetric=true)
    n_epochs = length(segments)

    # Output preallocation

    mask = falses(n_epochs)

    delta_power = zeros(Float64, n_epochs)
    beta_power  = zeros(Float64, n_epochs)

    # We compute PSDs and band powers once.
    for i in 1:n_epochs
        psd = PSD(
            segments[i],
            fs,
            fs*5;
            overlap = 0.5,
            window_function = hanning
        )

        delta_power[i] = total_band_power(psd, 0.6, 4.6)
        beta_power[i]  = total_band_power(psd, 40.0, 60.0)
    end

    # Sliding window evaluation
    for i in 1:n_epochs

        start_idx, end_idx = epoch_window_bounds(
            n_epochs,
            i,
            window_length
        )

        # Compute local averages excluding center epoch

        n = 0
        delta_sum = 0.0
        beta_sum  = 0.0

        # We know start_idx, end_idx are valid. We avoid 
        # safety checks for optimization.
        @inbounds for j in start_idx:end_idx
            if j != i
                delta_sum += delta_power[j]
                beta_sum  += beta_power[j]
                n += 1
            end
        end

        # Safety check
        if n == 0
            continue
        end

        delta_local = delta_sum / n
        beta_local  = beta_sum  / n

        # Ratio test
        if (delta_power[i] / delta_local > delta_threshold) ||
           (beta_power[i]  / beta_local  > beta_threshold)

            mask[i] = true
        end
    end

    return mask
end

"""
    hjorth_artifacts(signal::TimeSeries; z_thresh=3.0, k=1, mask=nothing, verbose=true) -> BitVector

Detect artifact epochs based on Hjorth parameters using an iterative
global outlier detection method, optionally restricted to a subset of epochs.

The signal is segmented into 30-second epochs. Hjorth parameters (Activity,
Mobility, Complexity) are computed for each epoch.

The outlier detection runs for `k` iterations. In each iteration:
1. Mean and Standard Deviation are calculated using only epochs that are:
   - Included in the input `mask` (if provided).
   - NOT marked as artifacts in previous rounds.
2. Epochs deviating by more than `z_thresh` from these new statistics are
   added to the artifact mask.

Segmentation is done symmetrically, so if the signal has an ending epoch of <30
seconds, it will will be dropped. This may result in the output mask having a
length of `num_epochs(signal) - 1` instead of `num_epochs(signal)`.

# Arguments
- `signal::TimeSeries`: Input time series.

# Keyword Arguments
- `z_thresh::Real=3.0`: Z-score threshold for identifying outliers.
- `k::Int=1`: Number of iterations for outlier detection.
- `mask::Union{BitVector, Nothing}=nothing`: Optional mask. If provided,
  artifact detection is performed ONLY on epochs where `mask` is `true`.
  Epochs where `mask` is `false` are ignored and will never be marked as artifacts.
- `verbose::Bool=false`: If `true`, prints detailed metrics about the detection process.

# Returns
- `Vector{Bool}`: Boolean mask where `true` indicates an epoch classified as artifact.
"""
function hjorth_artifacts(
    signal::TimeSeries;
    z_thresh::Real = 3.0,
    k::Int = 1,
    mask::Union{BitVector, Nothing} = nothing,
    verbose::Bool = true
)
    # 1. Setup
    fs = signal.fs
    epoch_len = fs * 30
    segments = segment(signal.x, epoch_len; symmetric=true)
    n_epochs = length(segments)

    # Validate input mask length if provided
    if !isnothing(mask) && length(mask) != n_epochs
        throw(DimensionMismatch("Input mask length ($(length(mask))) does not match number of epochs ($n_epochs)"))
    end

    # If no mask provided, treat all epochs as valid candidates
    candidate_mask = isnothing(mask) ? trues(n_epochs) : copy(mask)
    n_analyzed = count(candidate_mask)

    if verbose
        println("\n--- Hjorth Artifact Detection (Iterative) ---")
        println("Total epochs: $n_epochs")
        println("Analyzed epochs (mask): $n_analyzed")
        println("Params: z_thresh=$z_thresh, k=$k")
        println("---------------------------------------------")
    end

    activity   = Vector{Float64}(undef, n_epochs)
    mobility   = Vector{Float64}(undef, n_epochs)
    complexity = Vector{Float64}(undef, n_epochs)

    # 2. Parameter Calculation
    @inbounds for i in 1:n_epochs
        h = hjorth_parameters(segments[i])
        activity[i]   = h.activity
        mobility[i]   = h.mobility
        complexity[i] = h.complexity
    end

    # 3. Iterative Outlier Detection
    artifacts_mask = falses(n_epochs)

    for iter in 1:k
        # Identify epochs to use for statistics
        stats_indices = candidate_mask .& (.!artifacts_mask)
        
        # Safety check: Need at least 2 points
        if count(stats_indices) < 2
            if verbose
                println("Iteration $iter: Stopped (Insufficient data for statistics).")
            end
            break
        end

        # Calculate statistics using valid epochs
        @views begin
            act_clean = activity[stats_indices]
            mob_clean = mobility[stats_indices]
            com_clean = complexity[stats_indices]
        end

        μ_act, σ_act = mean(act_clean), std(act_clean)
        μ_mob, σ_mob = mean(mob_clean), std(mob_clean)
        μ_com, σ_com = mean(com_clean), std(com_clean)

        new_artifacts_count = 0
        found_new_artifact = false
        
        @inbounds for i in 1:n_epochs
            # Only check epochs that are candidates and not already marked
            if !candidate_mask[i] || artifacts_mask[i]
                continue
            end

            act_out = abs(activity[i] - μ_act) > z_thresh * σ_act
            mob_out = abs(mobility[i] - μ_mob) > z_thresh * σ_mob
            com_out = abs(complexity[i] - μ_com) > z_thresh * σ_com

            if act_out || mob_out || com_out
                artifacts_mask[i] = true
                found_new_artifact = true
                new_artifacts_count += 1
            end
        end

        if verbose
            println("Iteration $iter: flagged $new_artifacts_count new artifacts.")
        end

        if !found_new_artifact
            if verbose
                println("Converged at iteration $iter (no new artifacts found).")
            end
            break
        end
    end

    if verbose
        total_flagged = count(artifacts_mask)
        pct = n_analyzed > 0 ? round(100 * total_flagged / n_analyzed; digits=2) : 0.0
        
        println("---------------------------------------------")
        println("Final Results:")
        println("  Total Flagged: $total_flagged / $n_analyzed")
        println("  Percentage:    $pct%")
        println("---------------------------------------------\n")
    end

    return artifacts_mask
end
