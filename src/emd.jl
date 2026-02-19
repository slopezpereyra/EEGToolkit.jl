
"""
    find_extrema(x::AbstractVector)

Identifies the indices of local maxima and minima in a signal `x`. Explicitly
includes the endpoints (1 and N) to define the envelope over the full domain.
"""
function find_extrema(x::AbstractVector{T}) where T <: Real
    n = length(x)
    # Start with the first index to anchor the envelope
    max_inds = Int[1]
    min_inds = Int[1]
    
    # Scan for local extrema (strict inequality for simple peak detection)
    for i in 2:n-1
        if x[i-1] < x[i] > x[i+1]
            push!(max_inds, i)
        elseif x[i-1] > x[i] < x[i+1]
            push!(min_inds, i)
        end
    end
    
    # End with the last index to anchor the envelope
    push!(max_inds, n)
    push!(min_inds, n)
    
    return max_inds, min_inds
end

"""
    empirical_mode_decomposition(signal::AbstractVector; 
                                 max_imfs::Int=10, 
                                 sift_thresh::Float64=0.2, 
                                 max_sifts::Int=20)

Performs EMD using cubic splines.

# Arguments
- `signal`: Input vector.
- `max_imfs`: Safety limit for number of IMFs.
- `sift_thresh`: Standard Deviation threshold for stopping the sifting process.
- `max_sifts`: Max iterations per IMF.

"""
function emd(signal::AbstractVector{T}; 
                                      max_imfs::Int=10, 
                                      sift_thresh::Float64=0.2, 
                                      max_sifts::Int=20) where T <: AbstractFloat
    
    residue = copy(signal)
    imfs = Vector{Vector{T}}()
    n = length(signal)
    t_axis = 1:n  # The time axis for evaluation
    
    for imf_count in 1:max_imfs
        
        h = copy(residue)
        
        # Sifting Loop
        sift_sd = Inf
        sift_iter = 0
        
        while sift_sd > sift_thresh && sift_iter < max_sifts
            sift_iter += 1
            
            # 1. Find Extrema
            max_inds, min_inds = find_extrema(h)
            
            # Safety break: need at least 3 points (including edges) to define a cubic spline
            if length(max_inds) < 3 || length(min_inds) < 3
                break 
            end
            
            # 2. Interpolate Envelopes (Irregular Grid) DataInterpolations.jl
            # signature: CubicSpline(values, knots) We broadcast the result over
            # t_axis to get the envelope arrays.
            
            # Note: We convert indices to Float64 for type stability in math
            # operations, though DataInterpolations handles Int knots well.
            itp_up = CubicSpline(h[max_inds], Float64.(max_inds))
            itp_low = CubicSpline(h[min_inds], Float64.(min_inds))
            
            upper_env = itp_up.(t_axis)
            lower_env = itp_low.(t_axis)
            
            # 3. Mean Envelope & Subtraction
            mean_env = (upper_env .+ lower_env) ./ 2.0
            prev_h = copy(h)
            h = h .- mean_env
            
            # 4. Stop Criterion (SD method)
            eps_safe = 1e-10
            num = sum((prev_h .- h).^2)
            den = sum(prev_h.^2) + eps_safe
            sift_sd = num / den
        end
        
        # Check if this IMF is valid (contains extrema)
        max_inds, min_inds = find_extrema(h)
        if length(max_inds) < 3 || length(min_inds) < 3
            break
        end

        push!(imfs, h)
        residue = residue .- h
    end
    
    push!(imfs, residue)
    return imfs
end

"""
    hilbert_transform(imfs::Vector{Vector{T}}, fs::Real)

Applies the Hilbert transform to each IMF to extract instantaneous 
amplitude and frequency.

# Arguments
- `imfs`: Vector of IMFs (output from your `emd` function).
- `fs`: Sampling frequency in Hz.

# Returns
- `inst_amp`: Vector of vectors containing instantaneous amplitudes.
- `inst_freq`: Vector of vectors containing instantaneous frequencies (Hz).
"""
function hht(imfs::AbstractVector{<:AbstractVector{T}}, fs::Real) where T <: AbstractFloat
    
    n_imfs = length(imfs)
    n_samples = length(imfs[1])
    dt = 1.0 / fs
    
    # Pre-allocate output containers
    inst_amp = Vector{Vector{T}}(undef, n_imfs)
    inst_freq = Vector{Vector{T}}(undef, n_imfs)
    
    for (i, imf) in enumerate(imfs)
        # 1. Compute Analytic Signal using DSP.hilbert
        # z(t) = imf(t) + i * H[imf(t)]
        analytic_signal = hilbert(imf)
        
        # 2. Instantaneous Amplitude
        # A(t) = |z(t)|
        inst_amp[i] = abs.(analytic_signal)
        
        # 3. Instantaneous Phase
        # theta(t) = angle(z(t))
        # We must 'unwrap' the phase to avoid discontinuities at +/- pi
        theta = unwrap(angle.(analytic_signal))
        
        # 4. Instantaneous Frequency
        # f(t) = (1/2pi) * d(theta)/dt
        # We use simple finite differences for the derivative.
        # This reduces the array length by 1, so we pad the end to match size.
        d_theta = diff(theta)
        freq = d_theta ./ (2π * dt)
        
        # Pad the last point to maintain vector length (repeat last freq)
        push!(freq, freq[end])
        
        inst_freq[i] = freq
    end
    
    return inst_amp, inst_freq
end

# Plotting functions
#

"""
    plot_imfs(signal::AbstractVector, imfs::AbstractVector{<:AbstractVector})

Produces a plot of IMFs derived from the `emd` function.
- Top subplot: Original Signal.
- Middle subplots: Intrinsic Mode Functions (IMFs).
- Bottom subplot: Residual (Trend).

Returns the plot object for further manipulation if needed.
"""
function plot_imfs(signal::AbstractVector, imfs::AbstractVector{<:AbstractVector})
    n_imfs = length(imfs)
    # Total plots = 1 (Original) + n_imfs (IMFs + Residual)
    total_plots = n_imfs + 1
    
    # Define a vertical layout
    l = @layout [grid(total_plots, 1)]
    
    # Initialize plot with dynamic height (150px per component)
    p = plot(layout=l, size=(800, 150 * total_plots), legend=false)
    
    # 1. Plot Original Signal
    plot!(p[1], signal, title="Original Signal", linewidth=2, color=:black)
    
    # 2. Plot IMFs and Residual
    for (i, component) in enumerate(imfs)
        subplot_idx = i + 1
        
        # Check if it is the last component (Residual) or a standard IMF
        is_residual = (i == n_imfs)
        
        label = is_residual ? "Residual (Trend)" : "IMF $i"
        line_color = is_residual ? :red : :blue
        
        plot!(p[subplot_idx], component, 
              title=label, 
              linewidth=1.5, 
              color=line_color)
    end
    return p
end

"""
    plot_hilbert_spectrum(inst_amp, inst_freq, t_range; freq_lims=(0, 60))

Plots the Hilbert Spectrum with optional frequency constraints.

# Arguments
- `freq_lims`: A tuple (min_freq, max_freq) to filter visibility. 
               Default is (0, 60) which is standard for EEG.
               Pass `nothing` to see everything.
"""
function plot_hilbert_spectrum(inst_amp, inst_freq, t_range; freq_lims=(0, 60))
    
    # 1. Setup the plot limits
    # If limits are provided, use them. Otherwise, let Plots.jl auto-scale.
    if isnothing(freq_lims)
        p = plot(title="Hilbert Spectrum", xlabel="Time (s)", ylabel="Frequency (Hz)", legend=false)
    else
        p = plot(title="Hilbert Spectrum", xlabel="Time (s)", ylabel="Frequency (Hz)", 
                 legend=false, ylim=freq_lims)
    end
    
    n_imfs = length(inst_amp)
    
    # Loop up to n_imfs - 1 (skip residual)
    for i in 1:(n_imfs - 1)
        freq = inst_freq[i]
        amp = inst_amp[i]
        
        # 2. Create a Mask for Valid Frequencies
        # We allow a small buffer (< 0) just to catch numerical zeros, 
        # but strictly filter the upper bound if provided.
        if isnothing(freq_lims)
            mask = freq .> 0
        else
            (min_f, max_f) = freq_lims
            mask = (freq .>= min_f) .& (freq .<= max_f)
        end
        
        # Only plot if we have data points in this range
        if count(mask) > 0
            scatter!(p, 
                t_range[mask], 
                freq[mask], 
                zcolor=amp[mask],
                marker=(:circle, 2, stroke(0)),
                c=:viridis,
                colorbar_title="Amplitude"
            )
        end
    end
    
    display(p)
end

"""
    plot_hilbert_heatmap(inst_amp, inst_freq, t_range; 
                         freq_lims=(0, 30), 
                         time_bins=500, 
                         freq_bins=100)

Converts the scattered HHT data into a 2D Histogram (Heatmap).
This is much cleaner for long EEG recordings.

# Arguments
- `time_bins`: Number of vertical slices (resolution in time).
- `freq_bins`: Number of horizontal slices (resolution in frequency).
"""
function plot_hilbert_heatmap(inst_amp, inst_freq, t_range; 
                                  freq_lims=(0, 30), 
                                  res_t=800, 
                                  res_f=120)
    
    # 1. Define Edges
    t_min, t_max = t_range[1], t_range[end]
    f_min, f_max = freq_lims
    
    # 2. Pre-allocate the grid (Frequency x Time)
    spectrogram_matrix = zeros(Float64, res_f, res_t)
    
    dt_bin = (t_max - t_min) / res_t
    df_bin = (f_max - f_min) / res_f
    
    # 3. Loop and Accumulate (No Append!)
    # We iterate through each IMF and add its contribution to the grid
    n_imfs = length(inst_amp)
    
    for i in 1:(n_imfs - 1) # Skip residual
        freqs = inst_freq[i]
        amps = inst_amp[i]
        
        # Only iterate points within bounds (simple masking)
        for j in eachindex(freqs)
            f = freqs[j]
            t = t_range[j]
            a = amps[j]
            
            # Boundary check
            if f >= f_min && f < f_max && t >= t_min && t < t_max
                # Map continuous value to bin index
                t_idx = floor(Int, (t - t_min) / dt_bin) + 1
                f_idx = floor(Int, (f - f_min) / df_bin) + 1
                
                # Accumulate Energy (Amplitude^2)
                spectrogram_matrix[f_idx, t_idx] += a^2
            end
        end
    end
    
    # 4. Plot Heatmap
    heatmap(
        range(t_min, t_max, length=res_t),
        range(f_min, f_max, length=res_f),
        spectrogram_matrix,
        title="Fast Hilbert Spectrogram",
        xlabel="Time (s)", ylabel="Frequency (Hz)",
        c=:viridis,
        clims=(0, quantile(vec(spectrogram_matrix), 0.98)) # Auto-contrast
    )
end
