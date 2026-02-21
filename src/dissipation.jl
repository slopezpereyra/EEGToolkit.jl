# dissipation.jl

"""
    sw_dissipation(S::Spectrogram; kwargs...) -> (amp, k, errors)

Quantifies the homeostatic decline of Slow Wave Activity (SWA) by fitting a
2-parameter exponential decay model to the relative delta power time course.
This replicates the methodology of normalizing SWA relative to total NREM
amplitude as described in Armitage et al. (2000), "Slow-wave activity in NREM
sleep: sex and age effects in depressed outpatients and healthy controls".

The dissipation is modeled using a non-linear least squares approach as: 

SWA(t) = amp * exp(-k * t)

where `k` is the dissipation constant representing the rate of homeostatic decay.

# Arguments
- `S::Spectrogram`: The spectrogram object containing spectral data and the time domain.

# Keyword Arguments
- `initial_epoch::Int=1`: The starting epoch (spectrogram row) index for the analysis.
- `final_epoch::Int=length(S.time)`: The ending epoch (spectrogram row) index for the analysis.
- `start_params::Union{Nothing, Vector{Float64}}=nothing`: Initial guesses for 
  `[amplitude, k]`. If `nothing`, the function uses automatic heuristics.

# Process
1. **Trend Extraction**: Extracts absolute power in the delta band (0.4–4.0 Hz).
2. **Normalization**: Expresses SWA amplitude as a percentage of the total SWA amplitude 
   across the selected NREM epochs.
3. **Time Synchronization**: Shifts the time domain so that the fit starts at time = 0.
4. **Non-linear Fitting**: Utilizes least squares fitting to optimize the 2 parameters.

# Returns
- `amp::Float64`: The fitted amplitude (predicted SWA % at time = 0).
- `k_value::Float64`: The dissipation constant (decay rate).
- `errors::Vector{Float64}`: The standard errors for the fitted parameters.
"""
function sw_dissipation(S::Spectrogram; 
                        initial_epoch::Int=1, 
                        final_epoch::Int=length(S.time),
                        start_params::Union{Nothing, Vector{Float64}}=nothing)

    # 1. Validation 
    if initial_epoch < 1 || final_epoch > length(S.time) || initial_epoch >= final_epoch
        error("Invalid epoch range. Ensure 1 <= initial_epoch < final_epoch <= length(S.time)")
    end

    # 2. Extract absolute delta power (0.4 - 4.0 Hz) per epoch
    delta_trend_abs = absolute_band_power(S, 0.4, 4.0)
    
    # 3. Data Slicing. Float64 conversion for numerical stability in fitting.
    sliced_trend = Float64.(delta_trend_abs[initial_epoch:final_epoch])
    tdata_raw = Float64.(S.time[initial_epoch:final_epoch])

    # 4. Normalization and Time Synchronization
    # Express as percentage of total NREM amplitude
    total_swa = sum(sliced_trend)
    ydata = (sliced_trend ./ total_swa) .* 100
    
    # Synchronize time to sleep onset (start at t = 0)
    tdata = tdata_raw .- tdata_raw[1]

    # 5. Automatic Guessing Logic for the 2-parameter model
    p0 = if isnothing(start_params)
        # SWA % typically starts above 100% at the beginning of the night
        amp_guess = ydata[1] 
        # Standard decay constant observed in healthy adult cohorts
        k_guess = 0.001 
        [amp_guess, k_guess]
    else
        start_params
    end

    # 6. Define the 2-Parameter Exponential Decay Model
    # p[1]: Amplitude (b), p[2]: Decay Constant (k)
    @. model(t, p) = p[1] * exp(-p[2] * t)

    # 7. Execute Fit using LsqFit
    fit = curve_fit(model, tdata, ydata, p0)
    
    # 8. Extract Results
    amp, k_value = coef(fit)
    errors = stderror(fit) 
    
    return amp, k_value, errors
end

"""
    plot_dissipation_fit(S::Spectrogram, amp::Real, k::Real) -> Plots.Plot

Visualizes the normalized Slow Wave Activity (SWA) data against the 
fitted 2-parameter exponential decay curve.

# Arguments
- `S::Spectrogram`: The spectrogram object containing spectral data and time domain.
- `amp::Real`: The fitted amplitude parameter.
- `k::Real`: The fitted dissipation constant (decay rate).

# Returns
- `Plots.Plot`: The generated plot object overlaying the empirical data and fit.
"""
function plot_dissipation_fit(S::Spectrogram, amp::Real, k::Real)
    # 1. Reconstruct the normalized data space
    delta_trend_abs = absolute_band_power(S, 0.4, 4.0)
    total_swa = sum(delta_trend_abs)
    ydata_normalized = (delta_trend_abs ./ total_swa) .* 100

    # 2. Synchronize time to start at 0
    tdata_raw = Float64.(S.time)
    tdata_shifted = tdata_raw .- tdata_raw[1]

    # 3. Calculate the curve points using the provided parameters
    y_fit = amp .* exp.(-k .* tdata_shifted)

    # 4. Generate the visualization
    p = scatter(tdata_shifted, ydata_normalized, 
            label="Normalized NREM SWA", 
            xlabel="Time since NREM onset (epochs)", 
            ylabel="SWA (% of total NREM)",
            title="Homeostatic SWA Dissipation",
            markerstrokewidth=0,
            alpha=0.5)

    plot!(p, tdata_shifted, y_fit, 
          label="Exponential Fit (k = $(round(k, digits=4)))", 
          linewidth=3, 
          color=:red)
          
    return p
end
