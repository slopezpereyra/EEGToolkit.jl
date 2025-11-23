using Test
using EEGToolkit
using Statistics

@testset "EEGToolkit Tests" begin

    @testset "Connectivity: Coherence & wPLI" begin
        # 1. Setup Synthetic Data
        # We create two signals at 10Hz with a known relationship
        fs = 100.0
        t = 0:1/fs:1.0  # 1 second
        n_samples = length(t)
        n_epochs = 20
        
        # Signal X: Pure Sine wave at 10Hz
        x = repeat(sin.(2 * pi * 10 * t), 1, n_epochs)
        
        # Signal Y: Signal X + 90 degree phase shift (cosine)
        # wPLI should be maximal (1.0) for 90 degree shift
        y = repeat(cos.(2 * pi * 10 * t), 1, n_epochs)

        # 2. Test wPLI
        # We expect high wPLI at 10Hz (index 11 roughly, depending on FFT freq bins)
        wpli, freqs = compute_wpli(x, y, fs)
        
        # Find index closest to 10Hz
        idx_10hz = argmin(abs.(freqs .- 10.0))
        
        @test freqs[idx_10hz] ≈ 10.0 atol=0.5
        @test wpli[idx_10hz] > 0.9  # Should be very close to 1.0
        @test length(wpli) == length(freqs)

        # 3. Test Coherence
        # Coherence should also be 1.0 because they are perfectly linearly correlated
        coh, freqs_c = compute_coherence(x, y, fs)
        
        @test coh[idx_10hz] > 0.9
        @test all(coh .>= 0.0) && all(coh .<= 1.0 + eps()) # Bounds check
    end

end
