
"""
A struct representing time series data.

# Fields
- `x::Vector{<:AbstractFloat}`: Time series data.
- `fs::Integer`: Sampling rate.
"""
struct TimeSeries 
  x::Vector{<:AbstractFloat}
  fs::Integer 

  function TimeSeries(x, fs)
    new(x, fs)
  end
end

"""Given an epoch number, a sampling rate, and an epoch length in seconds, 
returns all indexes which would correspond to said epoch in a signal."""
function epoch_to_indexes(n::Integer, fs::Integer, epoch_length::Integer)::Tuple{Integer, Integer}
  s = ((n - 1) * fs * epoch_length) + 1
  e = n * fs * epoch_length
  (s, e)
end

"""Given two epoch numbers, a sampling rate, and an epoch length in seconds, 
returns all indexes which would correspond to the closed interval between both 
epochs in a signal."""
function epochs_to_indexes(n::Integer, m::Integer, fs::Integer, epoch_length::Integer)::Tuple{Integer, Integer}
  s = ((n - 1) * fs * epoch_length) + 1
  e = m * fs * epoch_length
  (s, e)
end


"""
Returns the duration in seconds of a time series.
"""
function length_in_secs(ts::TimeSeries)
  length(ts.x)/ts.fs
end

"""
Returns the duration in minutes of a time series.
"""
function length_in_mins(ts::TimeSeries)
  length_in_secs(ts)/60
end

function num_epochs(ts::TimeSeries, epoch_length::Int)::Int
  total_duration = length(ts.x) / ts.fs
  floor(Int, total_duration / epoch_length)
end

"""
Returns the duration in hours of a time series.
"""
function length_in_hours(ts::TimeSeries)
  length_in_mins(ts)/60
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
function segment(v::Vector{T}, L::Int; overlap::Union{<:AbstractFloat,Integer}=0, symmetric=false)::Vector{Vector{T}} where {T}
  if L > length(v)
    throw(ArgumentError("Segment length L must be less than or equal to the length of the vector."))
  end
  if overlap < 0 || overlap >= 1.0
    throw(ArgumentError("Overlap fraction must be in the range [0, 1)."))
  end

  if L == length(v)
    @warn("In the `segment` function, the length `L` of each segment was set to the length of `v`. Returning a vector containing `v`; id est, [ v ]`.")
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
function segment(ts::TimeSeries, L::Int; kargs...)::Vector{Vector{<:AbstractFloat}}
  segment(ts.x, L; kargs...)
end

"""
`seconds_to_time(seconds::Union{AbstractFloat, Integer})`

Helper function: maps a time in seconds to a Time object.
"""
function seconds_to_time(seconds::Union{AbstractFloat, Integer})::Time
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
julia> gen_time_domain(500, 10, 11)
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

julia> gen_time_domain(500, 600, 6000)
2700000-element Vector{Time}:
00:10:00.002
00:10:00.004
⋮
01:39:59.998
01:40:00
```
"""
function gen_time_domain(fs::Integer, s::Union{Integer, AbstractFloat}, e::Union{Integer, AbstractFloat})::Vector{Time}
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
function gen_time_domain(signal::TimeSeries, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer}; epoch_length::Integer = 30)::Vector{Time}
  gen_time_domain(signal.fs, s*epoch_length, e*epoch_length)
end

"""
`gen_time_domain(signal::TimeSeries)`

Generates a vector of Time objects representing the time instances ``t_1, \\ldots, t_n`` in
a `TimeSeries` signal, starting at `00:00:init` where `init` is ``\\frac{1}{f_s}``.
"""
function gen_time_domain(signal::TimeSeries)::Vector{Time}
  gen_time_domain(signal.fs, 0, length(signal.x) / signal.fs)
end


"""
`epoch(signal::TimeSeries, n::Integer; epoch_length::Integer=30)`

Returns a vector `[x₁, …, xₖ]` with all values `xᵢ` corresponding to the `n`th epoch in the signal.
"""
function epoch(signal::TimeSeries, n::Integer; epoch_length::Integer = 30)::TimeSeries
  range = epoch_to_indexes(n, signal.fs, epoch_length)
  if  range[1] < 0 || range[2] > length(signal.x)
    msg = "The epoch provided to the `epoch` function does not exist in the 
      signal, because it exceeds its duration or is negative"
    throw(ArgumentError(msg))
  end
  y = signal.x[range[1]:range[2]]
  TimeSeries(y, signal.fs)
end

"""
`epoch(signal::TimeSeries, n::Integer, m::Integer)`

Returns a vector `[x₁, …, xₖ]` with all indexes corresponding to epochs `n, n+1, …, m` of the EEG.
The default sampling rate is used to compute the indexes.
"""
function epoch(signal::TimeSeries, n::Integer, m::Integer; epoch_length::Integer = 30)::TimeSeries
  if (n == m)
    return epoch(signal, n)
  end
  if (n > m)
    throw(ArgumentError("The second epoch should be greater than the first."))
  end
  range = epochs_to_indexes(n, m, signal.fs, epoch_length)
  if  range[1] < 0 || range[2] > length(signal.x)
    msg = "The epoch provided to the `epoch` function does not exist in the 
      signal, because it exceeds its duration or is negative"
    throw(ArgumentError(msg))
  end
  y = signal.x[range[1]:range[2]]
  TimeSeries(y, signal.fs)
end


"""
`plot_ts(ts::TimeSeries, s::Integer, e::Integer; norm=false, ylab="Amplitude (uV)")`

Plots `TimeSeries` from epoch `s` to epoch `e`. The series many be normalized.
"""
function plot_ts(ts::TimeSeries, s::Integer, e::Integer; norm=false, ylab="Amplitude (uV)")::Plots.Plot
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
function plot_ts(ts::TimeSeries, s::Integer; kargs...)::Plots.Plot
  plot_ts(ts, s, s; kargs...)
end

"""
`plot_ts(ts::TimeSeries; norm=false, ylab="Amplitude (uV)")`

Plots `TimeSeries`. The series may be normalized.
"""
function plot_ts(ts::TimeSeries; norm=false, ylab="Amplitude (uV)")::Plots.Plot
  t = gen_time_domain(ts)
  y = norm ? ts.x .- mean(ts.x) ./ std(ts.x) : ts.x
  p = plot(t, y, ylabel = ylab, xlabel = "Time");
  return(p)
end

"""
function artifact_reject(signal::TimeSeries, anom_dict::Dict{Int, Vector{Int}}; epoch_length::Integer=30, subepoch_length::Integer=5)::Vector{Vector{<:AbstractFloat}};

This function removes from a signal the sub-epochs which contain artifacts. 
It requires a `TimeSeries` and a `Dict{Int, Vector{Int}}`, hereby termed 
`anom_dict` (for anomaly dictionary). 

`anom_dict` is understood to be such that `anom_dict[i] = [n₁, …, nₖ]`
means the `i`th epoch has artifacts at sub-epochs `n₁, ..., nₖ`.

The return value is a segmented signal (`Vector{Vector<:AbstractFloat}}`), each 
of whose segments corresponds to an epoch with its artifcat-contaminated 
sub-epochs removed. In other words, if the `result`   holds   the   return   
value   of    this    function,    `result[i]`
contains   what   is   left   from    the    `i`th    epoch    after    removing
its   contaminated   sub-epochs.   It   is   possible   that   `result[i]`    is
empty, if all sub-epochs of epoch `i` contained artifacts.
"""
function artifact_reject(signal::TimeSeries, anom_dict::Dict{Int, Vector{Int}}; epoch_length::Integer=30, subepoch_length::Integer=5)::Vector{Vector{<:AbstractFloat}}

  if any(x -> x < 0 || x > num_epochs(signal, epoch_length), keys(anom_dict))
    msg = "The `anom_dict` argument contains keys greater than the number of epochs in the `signal` or 
      non-positive."
    throw( ArgumentError(msg)  )
  end

  T = typeof(signal.x[1])

  epochs = segment(signal, signal.fs * epoch_length; symmetric = true)
  windows = map(e -> segment(e, signal.fs * subepoch_length; symmetric=true), epochs)

  # Iterate through each ( epoch , contaminated sub-epochs ) pair
  for (epoch, subepochs) in anom_dict
    # Sort the indexes of contaminated sub-epochs to be used in `deleteat!`.
    sorted_contaminated_epochs = sort(subepochs)
    # Delete from the epoch the contaminated sub-epochs.
    deleteat!(windows[epoch], sorted_contaminated_epochs)
  end

  # Collect each epoch after artifact removal, return vector of clean epochs.
  [Vector{T}(vcat(window...)) for window in windows]
end


"""
`artifact_reject(signal::TimeSeries, anoms::Vector{Integer}; epoch_length::Integer=30)`

An anomaly vector ``\\vec{x} \\in \\mathbb{N}^{n}`` is a sorted vector
whose values are those epochs in an `TimeSeries` that contain 
anomalies or artifacts. This function segments the TimeSeries and filters out all 
epochs containing artifacts.
"""
function artifact_reject(signal::TimeSeries, anoms::Vector{Integer}; epoch_length::Integer=30)
  epochs = segment(signal, signal.fs * epoch_length; symmetric = true)
  deleteat!( epochs, anoms )
end
