

"""
A struct representing time series data.

# Fields
- `x::Vector{<:AbstractFloat}`: Time series data.
- `fs::Integer`: Sampling rate.
- `epoch_length::Integer`: Length in seconds understood to comprise an epoch.
"""
struct TimeSeries 
    x::Vector{<:AbstractFloat}
    fs::Integer 
    epoch_length::Integer
end

"""
`segment(v::Vector{T}, L::Int, overlap_frac::Union{Float64,Int}) where {T}`

Splits a vector `v` into segments of length `L` with an overlap `overlap_frac` expressed as a fraction of L. 
"""
function segment(v::Vector, L::Int, overlap_frac::Union{Float64,Int})
    if L > length(v)
        throw(ArgumentError("Segment length L must be less than or equal to the length of the vector."))
    end

    if overlap_frac < 0 || overlap_frac >= 1.0
        throw(ArgumentError("Overlap fraction must be in the range [0, 1)."))
    end

    D = L * overlap_frac
    M = Int(ceil((length(v) - L) / (L - D)))

    segments = Vector{Vector{T}}(undef, M)
    step = Int(floor((1 - overlap_frac) * L))  # Calculate step size

    for i in 1:M
        start_idx = 1 + (i - 1) * step
        end_idx = start_idx + L - 1

        # Ensure the last segment does not exceed the length of the vector
        if end_idx > length(v)
            break
        end

        segments[i] = v[start_idx:end_idx]
    end

    return segments
end

"""
`segment(ts::TimeSeries, L::Int, overlap_frac::Union{Float64, Int})`

Splits the vector `ts.x` of a `TimeSeries` `ts` into segments of length `L`
with an overlap `overlap_frac` expressed as a fraction of L. 
"""
function segment(ts::TimeSeries, L::Int, overlap_frac::Union{Float64, Int})
    segment(ts.x, L, overlap_frac)
end

"""
`gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})`

Given an TimeSeries signal, generates the time vector t₁, …, tₙ corresponding to 
the signal from time `s` to `e` in seconds.
"""
function gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})
    step = 1 / signal.fs
    [i for i in (s+step):step:e]
end

"""
`epoch(signal::TimeSeries, n::Integer)`

Returns a vector `[x₁, …, xₖ]` with all values `xᵢ` corresponding to the `n`th epoch in the signal.
"""
function epoch(signal::TimeSeries, n::Integer)
    y = signal.x[((n - 1) * signal.fs * signal.epoch_length) + 1:n * signal.fs * signal.epoch_length]
    TimeSeries(y, signal.fs, signal.epoch_length)
end

"""
`epoch(signal::TimeSeries, n::Integer, m::Integer)`

Returns a vector `[x₁, …, xₖ]` with all indexes corresponding to epochs `n, n+1, …, m` of the EEG.
The default sampling rate is used to compute the indexes.
"""
function epoch(signal::TimeSeries, n::Integer, m::Integer)
    if (n == m)
        return epoch(eeg, n)
    end
    if (n > m)
        throw(ArgumentError("The second epoch should be greater than the first."))
    end
    y = signal.x[((n - 1) * signal.fs * signal.epoch_length) + 1:m * signal.fs * signal.epoch_length]
    TimeSeries(y, signal.fs, signal.epoch_length)
end
