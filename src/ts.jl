using Dates

"""
A struct representing time series data.

# Fields
- `x::Vector{<:AbstractFloat}`: Time series data.
- `fs::Integer`: Sampling rate.
- `secs::AbstractFloat`: Duration in seconds.
- `mins::AbstractFloat`: Duration in minutes.
- `hours::AbstractFloat`: Duration in hours.
- `epoch_length::Integer`: Length in seconds understood to comprise an epoch (defaults to 30).
- `subepoch_length::Integer`: Length in seconds understood to comprise an epoch (defaults to 30).
"""
struct TimeSeries 
    x::Vector{<:AbstractFloat}
    fs::Integer 
    secs::AbstractFloat
    mins::AbstractFloat
    hours::AbstractFloat
    epoch_length::Integer
    subepoch_length::Integer

    function TimeSeries(x, fs; epoch_length=30, subepoch_length=5)
        if subepoch_length >= epoch_length 
            throw(ArgumentError("Subepoch length must be inferior to epoch length."))
        end
        secs = length(x)/fs
        mins = secs/60
        hours = mins/60
        new(x, fs, secs, mins, hours, epoch_length, subepoch_length)
    end
end

"""
`segment(v::Vector{T}, L::Int; overlap::Union{<:AbstractFloat,Integer}=0, symmetric=false) where {T}`

Splits a vector `v` into segments of length `L` with an overlap `overlap` expressed as a fraction of L. The `overlap` defaults to `0` (no overlap).
Returns a vector ``v`` of vectors - i.e. `Vector{Vector{T}}` - with ``\\vec{v_i}`` the ``i``th segment in the split.

The function always attempts to capture the whole vector, even if the final split is not of length L. For example, 

```julia
> x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
> segment(x, 5)
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 4, 5]
 [6, 7, 8, 9, 0]

 > segment(x, 7)
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 4, 5, 6, 7]
 [8, 9, 0]
```

Set `symmetric=true` to ensure that, if this occurs, the last split is dropped.

```julia
> segment(x, 3; symmetric=true)
3-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [4, 5, 6]
 [7, 8, 9]
```

If `L` is equal to the segment length, `segment` raises a warning
and returns a vector with only the original vector: `[v]`. 
The return value ensures type-safety but the warning is raised because 
splitting a vector over its length is potentially 
a programming mistake.
"""
function segment(v::Vector{T}, L::Int; overlap::Union{<:AbstractFloat,Integer}=0, symmetric=false) where {T}
    if L > length(v)
        throw(ArgumentError("Segment length L must be less than or equal to the length of the vector."))
    end
    if overlap < 0 || overlap >= 1.0
        throw(ArgumentError("Overlap fraction must be in the range [0, 1)."))
    end

    if L == length(v)
        @warn("In the `segment` function, the length `L` of each segment was set to the length of `v`. Returning `v`.")
        return [v]
    end


    # Step size
    step = Int(floor(L - L * overlap))  # Calculate step size
    
    # Calculate number of segments
    M = Int(ceil((length(v) - L) / step)) + 1
    segments = Vector{Vector{T}}(undef, M)

    for i in 1:M
        start_idx = 1 + (i - 1) * step
        end_idx = start_idx + L - 1

        # Ensure the last segment does not exceed the length of the vector
        if start_idx > length(v)
            break
        end
        if end_idx > length(v)
            end_idx = length(v)
        end

        segments[i] = v[start_idx:end_idx]
    end

    # Filter out any undefined segments
    res = filter(!isnothing, segments)
    if symmetric && length( last(res) ) != L
        pop!(res)
    end
    
    return(res)
end

"""
`segment(ts::TimeSeries, L::Int; kargs...)`

Wrapper to segment the vector `ts.x` in the time 
series `ts`.
"""
function segment(ts::TimeSeries, L::Int; kargs...)
    segment(ts.x, L; kargs...)
end

"""
`seconds_to_time(seconds::Union{AbstractFloat, Integer})`

Helper function: maps a time in seconds to a Time object.
"""
function seconds_to_time(seconds::Union{AbstractFloat, Integer})
    # Calculate the components
    hours = floor(Int, seconds / 3600)
    minutes = floor(Int, (seconds % 3600) / 60)
    sec = floor(Int, seconds % 60)
    milliseconds = floor(Int, (seconds - floor(seconds)) * 1000)
    
    # Create the Time object
    Time(hours, minutes, sec, milliseconds)
end

"""
`gen_time_domain(fs::Integer, s::Union{Integer, AbstractFloat}, e::Union{Integer, AbstractFloat})`

Generates a vector of Time objects representing the time instances ``t_1, \\ldots, t_n`` in
a signal with with a given sampling rate ``f_s`` from second ``s`` to second ``e``. For instance,
``f_s = 500``, ``s = 10, e = 11`` would map to ``10.002, 10.004, \\ldots, 10.998, 11``.

```julia 
julia> gen_time_domain2(500, 10, 11)
500-element Vector{Time}:
 00:00:10.002
 00:00:10.004
 00:00:10.006
 00:00:10.008
 00:00:10.01
 00:00:10.012
 ⋮
 00:00:10.992
 00:00:10.994
 00:00:10.996
 00:00:10.998
 00:00:11

julia> gen_time_domain2(500, 600, 6000)
2700000-element Vector{Time}:
 00:10:00.002
 00:10:00.004
 ⋮
 01:39:59.998
 01:40:00
```
"""
function gen_time_domain(fs::Integer, s::Union{Integer, AbstractFloat}, e::Union{Integer, AbstractFloat})
    start_time_obj = seconds_to_time(s)
    step = 1 / fs
    L = length(collect(1:step:e-s+1))
    [start_time_obj + Millisecond(round(i * step * 1000)) for i in 1:L-1]
end

"""
`gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})`

Generates a vector of Time objects representing the time instances ``t_1, \\ldots, t_n`` in
a `TimeSeries` signal from epoch ``s`` to epoch ``e``. 
"""
function gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})
    gen_time_domain(signal.fs, s*signal.epoch_length, e*signal.epoch_length)
end

"""
`gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})`

Generates a vector of Time objects representing the time instances ``t_1, \\ldots, t_n`` in
a `TimeSeries` signal, starting at `00:00:init` where `init` is ``\\frac{1}{f_s}``.
"""
function gen_time_domain(signal::TimeSeries)
    gen_time_domain(signal.fs, 0, length(signal.x) / signal.fs)
end


"""
`epoch(signal::TimeSeries, n::Integer)`

Returns a vector `[x₁, …, xₖ]` with all values `xᵢ` corresponding to the `n`th epoch in the signal.
"""
function epoch(signal::TimeSeries, n::Integer)
    y = signal.x[((n - 1) * signal.fs * signal.epoch_length) + 1:n * signal.fs * signal.epoch_length]
    TimeSeries(y, signal.fs; epoch_length=signal.epoch_length, subepoch_length = signal.subepoch_length)
end

"""
`epoch(signal::TimeSeries, n::Integer, m::Integer)`

Returns a vector `[x₁, …, xₖ]` with all indexes corresponding to epochs `n, n+1, …, m` of the EEG.
The default sampling rate is used to compute the indexes.
"""
function epoch(signal::TimeSeries, n::Integer, m::Integer)
    if (n == m)
        return epoch(signal, n)
    end
    if (n > m)
        throw(ArgumentError("The second epoch should be greater than the first."))
    end
    y = signal.x[((n - 1) * signal.fs * signal.epoch_length) + 1:m * signal.fs * signal.epoch_length]
    TimeSeries(y, signal.fs; epoch_length=signal.epoch_length, subepoch_length = signal.subepoch_length)
end


"""
`plot_ts(ts::TimeSeries, s::Integer, e::Integer; norm=false, ylab="Amplitude (uV)") `

Plots `TimeSeries` from epoch `s` to epoch `e`. The series many be normalized.
"""
function plot_ts(ts::TimeSeries, s::Integer, e::Integer; norm=false, ylab="Amplitude (uV)")
    t = gen_time_domain(ts, s, e+1)
    ts = epoch(ts, s, e)
    y = norm ? ts.x .- mean(ts.x) ./ std(ts.x) : ts.x
    p = plot(t, y, ylabel = ylab, xlabel = "Time");
    return(p)
end

"""
`plot_ts(ts::TimeSeries, s::Integer; kargs...)`

Plots `TimeSeries` at epoch `s`.
"""
function plot_ts(ts::TimeSeries, s::Integer; kargs...)
    plot_ts(ts, s, s; kargs...)
end

"""
`plot_ts(ts::TimeSeries; norm=false, ylab="Amplitude (uV)")`

Plots `TimeSeries`. The series may be normalized.
"""
function plot_ts(ts::TimeSeries; norm=false, ylab="Amplitude (uV)")
    t = gen_time_domain(ts)
    y = norm ? ts.x .- mean(ts.x) ./ std(ts.x) : ts.x
    p = plot(t, y, ylabel = ylab, xlabel = "Time");
    return(p)
end


