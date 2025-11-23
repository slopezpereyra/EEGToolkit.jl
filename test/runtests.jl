using Pkg; Pkg.activate(".")

using EEGToolkit 
using Plots
using Statistics
using ContinuousWavelets
using LinearAlgebra
using Wavelets, FFTW
using DSP
using ImageFiltering

eeg = EEG("edfs/37bl.edf")
channel = get_channel(eeg, "EEG C4-A1")
signal = epoch(channel, 300).x
fs = 500.0
signal

# 1. FILTER THE SIGNAL (Crucial Step!)
# Real brain waves are mostly 0.5Hz - 30Hz. 
# Everything else is often noise (muscle, electrical hum) that creates that "static".
responsetype = Bandpass(0.5, 30; fs=fs)
designmethod = Butterworth(4)
clean_signal = filt(digitalfilter(responsetype, designmethod), signal)
# 1. Setup

# 1. Setup (Your signal)
# signal = ... (Already defined in your workspace)
n = length(clean_signal) 
fs = 500.0  # From your previous output

# 2. Define Wavelet (Exact syntax from example)
# We use Morlet(2π) because it handles EEG rhythms better than π
# Increase 'voices' (density of scales)
# default is usually 8-10. 32 gives a much smoother Y-axis.
c = wavelet(Morlet(2π), averagingType=NoAve(), β=2, voices=128)
# 3. Compute CWT
res = cwt(clean_signal, c)

# Parameters
epoch_idx = 300
epoch_len = 30.0 # seconds

# Calculate the start time of this specific epoch
start_time = (epoch_idx - 1) * epoch_len 

# Create the time vector shifted by start_time
# (0:(n-1)) ./ fs gives 0 to 30s
# Adding start_time shifts it to 8970s to 9000s
time_axis = start_time .+ (0:(n-1)) ./ fs

# 1. Compute Power
P = abs.(res).^2

# 2. Convert to Decibels
# We add eps() to avoid log(0)
p_db = 10 * log10.(P .+ eps())

# 3. CRITICAL STEP: Baseline Correction
# Subtract the median power of the whole signal to center the plot at 0 dB
p_db_centered = p_db .- median(p_db)

p_smooth = imfilter(p_db_centered, Kernel.gaussian((0.5, 2.0)))

# 4. Plot
heatmap(time_axis, 
    1:size(res, 2),
    p_smooth',       # Plot the CENTERED data
    xlabel = "Time (seconds)",
    ylabel = "Frequency Index", 
    title = "EEG CWT (dB relative to Median)",
    c = :jet,             
    clims = (-20, 20),     # Now these limits make sense!,
    interpolation = :cubic
)
