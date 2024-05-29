include("psd.jl")

function sigma_index(x::Vector{<:AbstractFloat}, fs::Integer)
    segs = overlaps(x, fs, 0)
    amps = map(x -> AmplitudeSpectrum(x, fs, 512), segs)

    # Helper function `f`: Applys function `g` to the amplitude spectrum in a given 
    # frequency band.
    function f(amp::AmplitudeSpectrum, band_lower::Number, band_upper::Number, g::Function)
        freq_band = findall(x -> x >= band_lower && x <= band_upper, amp.freq)
        freq_band_power = amp.spectrum[freq_band]
        g(freq_band_power)
    end
    
    # Helper function. Computes the σ-index of an AmplitudeSpectrum.
    function compute_index(amp::AmplitudeSpectrum)
        spindle_max_power =  f(amp, 10.5, 16, maximum)
        alpha_rejection_thresh = f(amp, 7.5, 10, maximum)
        if alpha_rejection_thresh > spindle_max_power 
            return 0
        end
        low_freq_mean_power = f(amp, 4, 10, y -> sum(y) / length(y))
        high_freq_mean_power = f(amp, 20, 40, y -> sum(y) / length(y))
        return 2 * spindle_max_power / ( low_freq_mean_power + high_freq_mean_power )
    end

    Σ = map(compute_index, amps)
    return(Σ)
end

function relative_spindle_power(x::Vector{<:AbstractFloat}, fs::Integer)
    segs = overlaps(x, fs, 0)
    amps = map(x -> AmplitudeSpectrum(x, fs, 512), segs)
    
    function f(amp::AmplitudeSpectrum, band_lower::Number, band_upper::Number, g::Function)
        freq_band = findall(x -> x >= band_lower && x <= band_upper, amp.freq)
        freq_band_power = amp.spectrum[freq_band]
        g(freq_band_power)
    end

    spindle_band_powers = map(x -> f(x, 11, 15, sum), amps)
    total_band_powers = map(x -> f(x, 0.5, 40, sum), amps)

    return spindle_band_powers ./ total_band_powers
end


