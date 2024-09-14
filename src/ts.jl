
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

function length_in_secs(ts::TimeSeries)
    length(ts.x)/ts.fs
end

function length_in_mins(ts::TimeSeries)
    length_in_secs(ts)/60
end

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


