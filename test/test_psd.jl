using DSP
using FFTW
using Statistics

# Tests for PSD structure
@testset "PSD Tests" begin
    # Sample data for testing
    fs = 100  # Sampling rate
    signal = sin.(2 * π * 5 * (0:1/fs:1)) + 0.5 * sin.(2 * π * 10 * (0:1/fs:1))  # 5Hz and 10Hz signals combined

    # Test direct PSD constructor
    @testset "Direct PSD Constructor" begin
        psd = PSD(signal, fs)
        @test length(psd.freq) == div(length(signal), 2) + 1
        @test psd.freq[1] == 0
        @test psd.freq[end] <= fs / 2
        @test psd.spectrum[1] >= 0
    end

    # Test PSD constructor with padding
    @testset "PSD with Padding" begin
        pad_length = next_power_of_two(length(signal))
        psd_padded = PSD(signal, fs, pad=pad_length)
        @test length(psd_padded.freq) == div(pad_length, 2) + 1
        @test psd_padded.spectrum[1] >= 0
    end

    # Test PSD constructor with segments (Welch method)
    @testset "Welch's Method" begin
        seg_length = 50
        psd_welch = PSD(signal, fs, seg_length)
        @test length(psd_welch.freq) == div(seg_length, 2) + 1
    end

end

# Tests for Spectrogram structure
@testset "Spectrogram Tests" begin
    # Sample data for testing
    fs = 100  # Sampling rate
    signal = sin.(2 * π * 5 * (0:1/fs:2)) + 0.5 * sin.(2 * π * 10 * (0:1/fs:2))  # 5Hz and 10Hz signals combined

    # Define a PSD function to pass to the spectrogram constructor
    function psd_function(seg)
        PSD(seg, fs)
    end

    # Test spectrogram constructor
    @testset "Spectrogram Constructor" begin
        window_length = 100
        spec = Spectrogram(signal, window_length, psd_function)
        @test size(spec.spectrums, 1) == div(length(signal), window_length)
        @test size(spec.spectrums, 2) == div(window_length, 2) + 1
    end

    # Test frequency band extraction from spectrogram
end

# Test helper functions
@testset "Helper Functions" begin
    # Test next_power_of_two
    @testset "Next Power of Two" begin
        @test next_power_of_two(5) == 8
        @test next_power_of_two(16) == 16
    end

    # Test zero_pad
    @testset "Zero Padding" begin
        v = [1.0, 2.0, 3.0]
        padded = zero_pad(v, 5)
        @test length(padded) == 5
        @test padded[4] == 0.0
    end
end
