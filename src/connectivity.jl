
"""
    compute_wpli(x, y, sampling_rate)

Computes the Weighted Phase Lag Index (wPLI) between two signals.

# Arguments
- `x::AbstractMatrix`: Signal X with shape (n_samples, n_epochs).
- `y::AbstractMatrix`: Signal Y with shape (n_samples, n_epochs).
- `sampling_rate::Number`: The sampling frequency in Hz.

# Returns
- `wpli::Vector`: The wPLI value for each frequency bin.
- `freqs::Vector`: The corresponding frequency labels (0 to Nyquist).
"""
function compute_wpli(x::AbstractMatrix, y::AbstractMatrix, sampling_rate::Number)
    n_samples, n_epochs = size(x)

    # 1. Compute the FFT along the time dimension (dim 1)
    # Hanning window applied to reduce spectral leakage
    win = hanning(n_samples)
    x_windowed = x .* win
    y_windowed = y .* win
    
    # FFT along dimension 1 (time)
    X_fft = fft(x_windowed, 1)
    Y_fft = fft(y_windowed, 1)

    # 2. Compute the Cross-Spectrum: Sxy = X * conj(Y)
    Sxy = X_fft .* conj.(Y_fft)

    # 3. Extract the Imaginary component: Im(Sxy)
    ImSxy = imag.(Sxy)

    # 4. Compute the Numerator: | E[ Im(Sxy) ] |
    # We take the mean across epochs (dims=2) to estimate expected value
    numerator = abs.(mean(ImSxy, dims=2))

    # 5. Compute the denominator: E[ | Im(Sxy) | ]
    denominator = mean(abs.(ImSxy), dims=2)

    # 6. Calculate wPLI
    # We add eps() to the denominator to prevent DivisionByZero errors
    wpli_spectrum = numerator ./ (denominator .+ eps())

    # --- Post-Processing for Output ---
    
    # We only care about the first half of frequencies (Positive Frequencies)
    # because the inputs are real signals.
    n_freqs = div(n_samples, 2) + 1
    
    # Consistent manual calculation (matches psd.jl style)
    freqs = rfftfreq(n_samples, sampling_rate)
    
    return vec(wpli_spectrum[1:n_freqs]), freqs
end

"""
    compute_coherence(x, y, sampling_rate)

Computes the Magnitude-Squared Coherence between two signals.

# Arguments
- `x::AbstractMatrix`: Signal X with shape (n_samples, n_epochs).
- `y::AbstractMatrix`: Signal Y with shape (n_samples, n_epochs).
- `sampling_rate::Number`: The sampling frequency in Hz.

# Returns
- `coherence::Vector`: The coherence value (0 to 1) for each frequency bin.
- `freqs::Vector`: The corresponding frequency labels (0 to Nyquist).
"""
function compute_coherence(x::AbstractMatrix, y::AbstractMatrix, sampling_rate::Number)
    n_samples, n_epochs = size(x)

    # 1. Compute the FFT along the time dimension (dim 1)
    # Hanning window applied to reduce spectral leakage
    win = hanning(n_samples)
    x_windowed = x .* win
    y_windowed = y .* win
    
    # FFT along dimension 1 (time)
    X_fft = fft(x_windowed, 1)
    Y_fft = fft(y_windowed, 1)

    # 2. Compute Cross-Spectrum and Auto-Spectra
    # Cross-Spectrum: X * conj(Y)
    Sxy = X_fft .* conj.(Y_fft)
    
    # Auto-Spectra (Power): |X|^2 and |Y|^2
    # Note: abs2(x) is more efficient than abs(x)^2
    Sxx = abs2.(X_fft)
    Syy = abs2.(Y_fft)

    # 3. Average across epochs (Expectation E[...])
    # It is crucial to average Sxy (complex) BEFORE taking the magnitude
    Pxy = mean(Sxy, dims=2) 
    Pxx = mean(Sxx, dims=2)
    Pyy = mean(Syy, dims=2)

    # 4. Compute Magnitude-Squared Coherence
    # Formula: |Pxy|^2 / (Pxx * Pyy)
    numerator = abs2.(Pxy)
    denominator = Pxx .* Pyy

    # Add eps() to prevent division by zero
    coh_spectrum = numerator ./ (denominator .+ eps())

    # --- Post-Processing for Output ---
    
    # Keep only positive frequencies (0 to Nyquist)
    n_freqs = div(n_samples, 2) + 1
    
    # Frequency vector
    freqs = rfftfreq(n_samples, sampling_rate)
    
    return vec(coh_spectrum[1:n_freqs]), freqs
end
