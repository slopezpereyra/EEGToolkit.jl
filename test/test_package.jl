using EEGToolkit
using Test


DESIRED_CHANNELS = 
    [
    "EEG C4-A1", 
    "EEG F3-A2",
    "EEG C3-A2",
    "EEG F4-A1",
    "EEG O2-A1",
    "EEG O1-A2"
    ]

SPECT1 = 
    [
    48.59907082063545,  
    42.62304733744895,
    18.755319639502822,
    25.384217345669775,
    22.523079690173702,
    21.01119324673764,
    27.014207551611694,
    26.001251503887918,
    17.554512690898243,
    22.613672931765475,
    ]

SPECT2 = 
    [
 23.742243551064213,   
 22.021681892103445,
 19.184606951843556,
 14.012680131490239,
 12.041065369086203,
 10.526440778196168,
 8.911171980700793, 
 6.495502765634059,
 5.540642734468156,
 4.893498970548876,
    ]




@testset "EEGToolkit.jl" begin
    ε = 0.001
    edf_file = "edfs/37bl.edf"
    eeg = EEG(edf_file)
   
    filter_channels!(eeg, p -> startswith(first(p), "EEG"))
    remaining_channels = get_Channels(eeg)
#    filter!(p -> startswith(first(p), "EEG"), eeg.signals)
    @test collect(keys(remaining_channels)) == DESIRED_CHANNELS

    signal = epoch(get_channel(eeg, "EEG C4-A1"), 10, 11) 
    fs = signal.fs

    seg = segment(signal, 100)
    @test length(seg) == 300 


    psd1 = PSD(seg[1], fs; dB=true)
    @test all(x -> x < ε, abs.(psd1.spectrum[1:10] .- SPECT1))
    psd2 = PSD(s.x, fs, fs*5; dB=true)
    @test all(x -> x < ε, abs.(psd2.spectrum[1:10] .- SPECT2))

    spec = Spectrogram(s.x, fs*30, x -> PSD(x, fs, fs*5); dB=true)
    @test abs(spec.spectrums[1,3] - 16.775304018633083) < ε
    @test abs(spec.spectrums[2,1200] - -51.80621574940826) < ε

end
