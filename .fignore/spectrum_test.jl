using Random
using Plots
include("src/ts.jl")
include("src/psd.jl")

# Set seed for reproducibility
Random.seed!(123)

# Parameters
sampling_rate = 500       # Sampling rate (samples per second)
duration = 5 - 1/sampling_rate   # Duration of the signal in seconds
frequencies = [50, 100, 150]   # Frequencies of the sine waves (in Hz)
amplitudes = [1, 2, 4]        # Amplitudes of the sine waves

# Generate time vector
t = 0:1/sampling_rate:duration
n = length(t)

# Initialize signal as zeros
signal = zeros(n)

# Create a signal as a mixture of sine waves
for i in 1:length(frequencies)
    frequency = frequencies[i]
    amplitude = amplitudes[i]
    
    # Generate the sine wave component
    sine_wave = amplitude * sin.(2 * π * frequency .* t)
    
    # Add the sine wave component to the signal
    signal .+= sine_wave
end


amp = AmplitudeSpectrum(signal, 500)
plot(amp.freq, amp.spectrum)

psd = PSD(signal, 500; dB=true)
plot_psd(psd; freq_lim=250)

wpsd = PSD(signal, 500, 500*2; normalization = 1, inner_normalization=1, dB=true)
plot_psd(wpsd; freq_lim=250)






