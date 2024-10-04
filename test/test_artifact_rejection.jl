
@testset "Artifact Rejection Tests" begin

    @testset "Basic functionality" begin

        # Test format: Generate a 90 element vector with sampling rate 10.
        # This gives 9 seconds.
        # Set epoch length to 3, subepoch length to 1.
        # This gives 3 epochs, each with 3 subepochs.
        signal_data = TimeSeries(collect(1.1:0.1:10.0), 10)  # Signal from 1.0 to 10.0, fs = 10 Hz
        
        anom_dict = Dict(1 => [1], 2 => [2])  # Remove first subepoch of first epoch, second subepoch of second epoch
        result = artifact_reject(signal_data, anom_dict, epoch_length=3, subepoch_length=1)

        x = collect(2.1:0.1:4.0)
        y₁ = collect(4.1:0.1:5.0); y₂ = collect(6.1:0.1:7.0); y = append!(y₁, y₂)
        z = collect(7.1:0.1:10.0)

        @test length(result) == 3  # Should still have 3 epochs
        @test result[1] == x  # First subepoch removed from the first epoch
        @test result[2] == y  # Second subepoch removed from the second epoch
        @test result[3] == z  # Second subepoch removed from the second epoch
    end

    @testset "All subepochs in an epoch are deleted" begin
        signal_data = TimeSeries(collect(1.0:0.1:10.0), 10)
        anom_dict = Dict(1 => [1, 2, 3])  # Remove all subepochs from first epoch

        result = artifact_reject(signal_data, anom_dict, epoch_length=3, subepoch_length=1)
        @test length(result) == 3
        @test result[1] == Float64[]  # Entire epoch 1 is deleted
    end

    @testset "No anomalies provided (normal case)" begin
        signal_data = TimeSeries(collect(1.0:0.1:10.0), 10)
        anom_dict = Dict()  # No anomalies
        @test_throws MethodError artifact_reject(signal_data, anom_dict; epoch_length=3, subepoch_length=1)
    end

    @testset "Whole signal gets deleted" begin
        signal_data = TimeSeries(collect(1.1:0.1:10.0), 10)
        anom_dict = Dict(1 => [1, 2, 3], 2 => [1, 2, 3], 3 => [1, 2, 3])  # Delete all epochs

        result = artifact_reject(signal_data, anom_dict, epoch_length=3, subepoch_length=1)
        @test all(isempty, result)  # All epochs should be empty
    end

    @testset "Signal shorter than epoch length" begin
        signal_data = TimeSeries(collect(1.0:0.1:2.5), 10)  # Signal is shorter than one full epoch
        anom_dict = Dict(1 => [1, 2, 3], 2 => [1, 2, 3], 3 => [1, 2, 3])  # Delete all epochs

        @test_throws ArgumentError artifact_reject(signal_data, anom_dict; epoch_length=3, subepoch_length=1)
    end

    @testset "Subepoch length greater than epoch length" begin
        signal_data = TimeSeries(collect(1.0:0.1:10.0), 10)
        anom_dict = Dict(1 => [1])  # Remove subepoch 1 from epoch 1

        @test_throws ArgumentError artifact_reject(signal_data, anom_dict, epoch_length=3, subepoch_length=4)  # subepoch > epoch
    end

    @testset "Anomaly refers to a non-existent epoch or subepoch" begin
        signal_data = TimeSeries(collect(1.0:0.1:10.0), 10)
        anom_dict = Dict(5 => [1])  # Refers to a non-existent epoch

        @test_throws ArgumentError artifact_reject(signal_data, anom_dict, epoch_length=3, subepoch_length=1)
    end

end
