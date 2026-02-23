 using Test
using EEGToolkit
using Statistics
using DSP


# Set a seed for reproducible simulated data if desired
# using Random; Random.seed!(42)

@testset "EEGToolkit.jl Test Suite" begin

    @testset "TimeSeries & Segmentation" begin
        # Simulate a 10-second signal at 100 Hz
        fs = 100
        x = randn(fs * 10)
        ts = TimeSeries(x, fs) #
        
        @test length_in_secs(ts) == 10.0 #
        @test num_epochs(ts, 2) == 5 #
        
        # Segmenting into 2-second windows (200 samples)
        segs = segment(ts, fs * 2) #
        @test length(segs) == 5
        @test length(segs[1]) == fs * 2
    end

    @testset "Filtering & Resampling" begin
        fs = 500
        t = 0:(1/fs):(5 - 1/fs)
        # Create a signal mixing 5 Hz (low) and 50 Hz (high)
        x = sin.(2 * π * 5 .* t) .+ 0.5 .* sin.(2 * π * 50 .* t)
        ts = TimeSeries(x, fs)

        # Lowpass filter at 20 Hz
        ts_low = apply_lowpass(ts, 20.0) #
        
        # The 50 Hz component should be heavily attenuated
        pow_original = sum(abs2, x)
        pow_filtered = sum(abs2, ts_low.x)
        @test pow_filtered < pow_original

        # Resampling to 100 Hz
        ts_resampled = EEGToolkit.resample(ts, 100.0)
        @test ts_resampled.fs == 100
        @test length(ts_resampled.x) == 500 # 5 seconds * 100 Hz
    end

    @testset "Staging & Masks" begin
        # 1-4 = NREM, 5 = REM, 6 = Wake
        stages = ["6", "1", "2", "2", "3", "3", "5", "6"]
        stg = Staging(stages) #
        
        @test length(stg) == 8
        
        # Stage mask for NREM
        nrem_mask = stage_mask(stg; include=:NREM) #
        @test sum(nrem_mask) == 5 # One "1", two "2"s, two "3"s
    end

    @testset "Spectral Analysis (PSD & Spectrogram)" begin
        fs = 100
        t = 0:(1/fs):(10 - 1/fs)
        # 10 Hz sine wave
        x = sin.(2 * π * 10 .* t)
        
        # PSD
        psd = PSD(x, fs) #
        @test length(psd.freq) == length(psd.spectrum)
        
        # Peak frequency should be exactly at or very close to 10 Hz
        peak_idx = argmax(psd.spectrum)
        @test isapprox(psd.freq[peak_idx], 10.0, atol=1.0)
        
        # Spectrogram: 2-second windows
        spec = Spectrogram(x, 2 * fs, sig -> PSD(sig, fs)) #
        @test size(spec.spectrums, 1) == 5 # 10 seconds / 2-sec windows
        @test size(spec.spectrums, 2) == length(spec.freq)
    end


    @testset "EMD & Hilbert-Huang Transform" begin
        fs = 100
        t = 0:(1/fs):(2 - 1/fs)
        x = sin.(2 * π * 4 .* t) .+ 0.5 .* sin.(2 * π * 15 .* t)
        
        # EMD Decomposition
        imfs = emd(x; max_imfs=3) #
        @test length(imfs) > 1
        @test length(imfs[1]) == length(x)
        
        # HHT
        inst_amp, inst_freq = hht(imfs, fs) #
        @test length(inst_amp) == length(imfs)
        @test length(inst_freq[1]) == length(x)
    end

    @testset "Connectivity Metrics" begin
        fs = 100
        n_samples = 200
        n_epochs = 5
        
        # Simulated multi-epoch data: X and Y are identical (perfect coherence)
        X = randn(n_samples, n_epochs)
        Y = copy(X)
        
        coh, freqs = compute_coherence(X, Y, fs) #
        wpli, _ = compute_wpli(X, Y, fs) #
        
        @test length(coh) == length(freqs)
        @test all(isapprox.(coh, 1.0, atol=1e-5)) # Identical signals = Coherence of 1.0
        
        # wPLI for identical signals without phase lag should theoretically be 0
        @test all(wpli .< 0.1) 
    end

    @testset "Slow Wave Detection" begin
        fs = 100
        t = 0:(1/fs):(4 - 1/fs)
        # Create a slow oscillation (~1 Hz) with a large amplitude to trigger detection
        x = -80.0 .* sin.(2 * π * 1.0 .* t) # Starts with a negative phase
        
        waves = detect_slow_waves(x, fs; amp_neg=40.0, amp_ptp=70.0) #
        
        # It should detect at least one slow wave meeting the criteria
        @test length(waves) >= 1
        if length(waves) >= 1
            @test waves[1].ptp_amp >= 70.0 #
            @test waves[1].duration > 0.5 #
        end
    end

end   


