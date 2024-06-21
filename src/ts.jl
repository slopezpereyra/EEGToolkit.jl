

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
`gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})`

Given an TimeSeries, generates the time vector t₁, …, tₙ corresponding to 
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
    TimeSeries(y, signal.fs; epoch_length=signal.epoch_length, subepoch_length = signal.subepoch_length)
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
    TimeSeries(y, signal.fs; epoch_length=signal.epoch_length, subepoch_length = signal.subepoch_length)
end

"""
`plot_eeg(eeg::EEG, s::Integer, e::Integer; channels::Vector{String}=[""], spacing::AbstractFloat=1.5)`

Plots EEG channels from epoch `s` to epoch `e`. Specific channels may be selected with the `channels` karg.
The `spacing` argument is an added factor in the normalization of the EEG signals - the vertical distance between
each signal in the plot grows proportionally to `spacing`.
"""
function plot_eeg(eeg::EEG, s::Integer, e::Integer; channels::Vector{String}=[""], spacing::AbstractFloat=1.5)

    if any(x -> x ∉ collect(keys(eeg.signals)), channels)
        throw(ArgumentError("Attempting to plot a non-existent channel. Revise the `channels` keyword argument."))
    end
   
    signal_dict = channels != [""] ? filter(p -> first(p) in channels, eeg.signals) : eeg.signals
    S = [] # Channel values
    C = collect(keys(signal_dict)) # Channel names
    i = 1
    L = length(C)
    for (key, value) in signal_dict
        signal = epoch2(value, s, e).x
        signal = i .+ (signal .- mean(signal)) ./ (spacing * L * std(signal))
        push!(S, signal) 
        i += 1
    end
    p = plot(ylabel = "Amplitude (uV)", xlabel = "Time (s)", yticks = (1:L, C));

    for (i, s) in enumerate(S)
        plot!(1:length(s), s, label = "")
    end
    return(p)
end

"""
`remove_channel!(eeg::EEG, channel::String)`

Removes a channel from the EEG.
"""
function remove_channel!(eeg::EEG, channel::String)
    delete!(eeg.signals, channel)
end

"""
`remove_channel!(eeg::EEG, channel::String)`

Removes a list of channels from the EEG.
"""
function remove_channel!(eeg::EEG, channels::Vector{String})
    for chan in channels 
        delete!(eeg.signals, chan)
    end
end
