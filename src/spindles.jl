
"""
`sigma_index(x::Vector{<:AbstractFloat}, fs::Integer)`

The ``\\sigma``-index algorithm [(Huupponen et al.,
2007)](https://pubmed.ncbi.nlm.nih.gov/17555950/) find abnormally high amplitude values in the spindle frequency
band. Per each 1 second window of the EEG, it computes

- the maximum amplitude in the spindle frequency band, which we call ``S_{max}``
- the average amplitude in the low alpha and theta frequencies, which we call
``\\alpha_{mean}, \\theta_{mean}``
- the maximum alpha amplitude ``\\alpha_{max}``

The ``\\sigma``-index of each window is defined to be zero if ``\\alpha_max > S_{max}``,
and otherwise

```math 
f(S_{max}, \\alpha_{mean}, \\phi_{mean}) = \\frac{2S_{max}}{ \\alpha_{mean} + \\theta_{mean} } 
```

Higher values are indicative of a higher spindle probability. The rejection
threshold recommended in the original paper is ``\\lambda = 4.5``.
"""
function sigma_index(x::Vector{<:AbstractFloat}, fs::Integer)
  segs = segment(x, fs)
  amps = map(x -> AmplitudeSpectrum(x, fs), segs)

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
    theta_mean = f(amp, 4, 8, y -> sum(y) / length(y))
    alpha_mean = f(amp, 8, 10, y -> sum(y) / length(y))
    return 2 * spindle_max_power / (alpha_mean + theta_mean)
  end

  Σ = map(compute_index, amps)
  return(Σ)
end

"""
`relative_spindle_power(x::Vector{<:AbstractFloat}, fs::Integer)`

The Relative Spindle Power (RSP) algorithm [(Devuyst et al., 2011)](https://pubmed.ncbi.nlm.nih.gov/22254656/) 
also detects abnormal values along the spindle frequency band. 
For every 1 second window, the amplitude spectrum ``S(t)`` is computed, 
and the RSP is defined as

```math 
RSP(t) = \\frac{\\int_{11}^{16} S(t, f) df}{\\int_{0.5}^{40} S(t, f) df}
```

This definition is more intelligible than the that of the ``\\sigma``-index, insofar
as it represents the ratio of the total power in the spindle band with respect
to the total power in the ``[0.5, 40]`` frequency range. It is evident
by definition that ``0 \\leq RSP \\leq 1``. Higher values are indicative of a higher spindle
probability (it should be clear that ``RSP`` is not a probability itself).
The rejection threshold recommended in the original paper is ``\\lambda = 0.22``.
"""
function relative_spindle_power(x::Vector{<:AbstractFloat}, fs::Integer)
  segs = segment(x, fs)
  amps = map(x -> AmplitudeSpectrum(x, fs, 512), segs)

  function f(amp::AmplitudeSpectrum, band_lower::Number, band_upper::Number, g::Function)
    freq_band = findall(x -> x >= band_lower && x <= band_upper, amp.freq)
    freq_band_power = amp.spectrum[freq_band]
    g(freq_band_power)
  end

  spindle_band_powers = map(x -> f(x, 11, 16, sum), amps)
  total_band_powers = map(x -> f(x, 0.5, 40, sum), amps)

  return spindle_band_powers ./ total_band_powers
end

# ----------------





























