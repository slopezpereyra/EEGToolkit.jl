# This file is merely a token to avoid automatic test execution on each push.
# Tests are run with a Makefile.

using EEGToolkit 
using Test 

include("test_ts.jl")
include("test_psd.jl")
include("test_package.jl")
include("test_artifact_rejection.jl")
include("test_nrem.jl")
