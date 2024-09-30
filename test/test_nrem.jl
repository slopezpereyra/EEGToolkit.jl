include("../src/nrem.jl")
using Test

# Helper function for running the tests
function run_nrem_tests()
    # Test 1: Simple case with stages 1, 2, 3, and REM
    @testset "Simple NREM detection" begin
        staging = ["1", "2", "3", "5", "5", "5", "2", "3", "5", "5", "5", "2", "3", "5", "5", "3", "3", "3", "5"]
        result = nrem(staging, 30, 3)
        @test isnothing(result)  # There should be 1 NREM period
        result = nrem(staging, 2, 3)
        @test result == [[2, 3], [7, 8], [12, 13, 16, 17, 18]]  
    end

    # Test 2: Case with interspersed wake and REM stages
    @testset "Interspersed wake and REM stages" begin
        staging = ["1", "2", "6", "6", "6", "1", "2", "5", "3", "5", "1", "2", "3", "6"]
        result = nrem(staging, 1, 3)
        @test length(result) == 2  # There should be 2 NREM periods
        @test result[1] == [2]  # First NREM period (before wake)
        @test result[2] == [7, 9, 12, 13]  # Second NREM period (after REM)
    end

    # Test 7: Invalid stages in the input
    @testset "Invalid input handling" begin
        # Using a non-valid stage in the staging array
        staging = ["1", "2", "X", "6", "1", "2", "5", "3", "5", "1", "2", "3", "6"]
        @test_throws ArgumentError nrem(staging, 1, 1)
    end
end

run_nrem_tests()
