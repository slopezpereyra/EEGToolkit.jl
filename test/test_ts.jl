include("../src/ts.jl")
using Dates
using Test
using Plots
using Statistics

function gen_ts()
    x = [0.1, 0.2, 0.3, 0.4, 0.5]  # Example time series data
    fs = 2  # Sampling rate
    ts = TimeSeries(x, fs)
    return ts
end

# Test the creation of a TimeSeries
@testset "TimeSeries tests" begin
    ts = gen_ts()
    
    @test ts.x == [0.1, 0.2, 0.3, 0.4, 0.5]  # Data is correct
    @test ts.fs == 2  # Sampling rate is correct
end

# Test length calculation functions
@testset "Length calculation tests" begin
    ts = gen_ts()
    
    @test length_in_secs(ts) == length(ts.x) / ts.fs  # Check length in seconds
    @test length_in_mins(ts) == length_in_secs(ts) / 60  # Check length in minutes
    @test length_in_hours(ts) == length_in_mins(ts) / 60  # Check length in hours
end

# Test segment function for vector input
@testset "Segment function tests" begin
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
    
    @test segment(x, 5) == [[1, 2, 3, 4, 5], [6, 7, 8, 9, 0]]  # Check simple split
    @test segment(x, 7) == [[1, 2, 3, 4, 5, 6, 7], [8, 9, 0]]  # Check uneven split
    @test segment(x, 3; symmetric=true) == [[1, 2, 3], [4, 5, 6], [7, 8, 9]]  # Check symmetric split
end

# Test segment function for TimeSeries input
@testset "TimeSeries segment function tests" begin
    ts = gen_ts()
    
    @test segment(ts, 2) == [[0.1, 0.2], [0.3, 0.4], [0.5]]  # Segmenting the TimeSeries data
end

# Test gen_time_domain for specific intervals
@testset "gen_time_domain tests" begin
    fs = 500
    start_time = 10
    end_time = 11
    
    time_domain = gen_time_domain(fs, start_time, end_time)
    
    @test length(time_domain) == fs  # The length of the generated time domain should match the sampling rate
    
    # Ensure the time domain contains Time objects and their correct sequence
    @test time_domain[1] == Time(0, 0, 10, 2)  # First time step is 10.002 seconds
    @test time_domain[end] == Time(0, 0, 11, 0)  # Last time step is exactly 11.000 seconds
end

# Test epoch extraction
@testset "Epoch tests" begin
    ts = gen_ts()
    
    # Test single epoch extraction
    epoch_ts = epoch(ts, 1; epoch_length=1)
    @test epoch_ts.x == [0.1, 0.2]  # First epoch is the first 30 seconds of data
    @test epoch_ts.fs == 2  # fₛ should remain unchanged
    @test_throws ArgumentError epoch(ts, 3; epoch_length=1)
    @test_throws ArgumentError epoch(ts, 0; epoch_length=1)
    
    # Test multiple epochs extraction
   
    ts = TimeSeries([float( i ) for i in range(1, 1000)], 10)  # Create a TimeSeries with 1000 data points, fs = 10Hz, 3 and something epohcs
    epoch_ts = epoch(ts, 1, 3)  # Extract epochs 1-3
    @test epoch_ts.x == collect(1:900)  # Should return 900 data points (30 seconds of data from each epoch)
    @test_throws ArgumentError epoch(ts, 1, 4)
end

# Test helper functions
@testset "Helper function tests" begin
    @test seconds_to_time(3661.123) == Time(1, 1, 1, 123)  # Convert seconds to Time object
end

@testset "plot_ts tests" begin
    # Create random data with 30
    ts = TimeSeries(rand(1200), 10)
    plot_ts(ts, 1, 2)
    # Test 1: Ensure the function runs without errors for a range
    @test plot_ts(ts, 1, 2) isa Plots.Plot

    # Test 2: Ensure the function runs without errors for a single epoch
    @test plot_ts(ts, 1, 2) isa Plots.Plot

    # Test 3: Ensure the function runs without errors with normalization
    @test plot_ts(ts, 1, 2; norm=true) isa Plots.Plot

    # Test 4: Ensure the function returns a plot for full series
    @test plot_ts(ts) isa Plots.Plot

    # Test 5: Check invalid parameters

end
