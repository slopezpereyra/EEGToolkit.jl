
"""
sigma_index(x::AbstractVector{<:AbstractFloat}, fs::Integer; mask::Union{Nothing, BitVector}=nothing, threshold::Real=4.5)

Computes the sigma-index (Huupponen et al., 2007) for 1-second segments of EEG data.
Returns a DataFrame containing only the segments where the sigma-index exceeds the `threshold`.

# Arguments
- `x`: The EEG signal.
- `fs`: Sampling frequency.
- `mask`: Optional BitVector corresponding to 30-sec epochs to include (true) or ignore (false).
- `threshold`: The minimum sigma-index to be considered a spindle (default: 4.5).

Reference: https://pubmed.ncbi.nlm.nih.gov/17555950/
"""
function sigma_index(ts::TimeSeries; mask::Union{Nothing, BitVector}=nothing, threshold::Real=4.5)
    x = ts.x
    fs = ts.fs
    segs = segment(x, fs; symmetric=true)

    if mask !== nothing
        expected_epochs = num_epochs(ts, 30)
        if length(mask) != expected_epochs
            throw(ArgumentError("Mask length ($(length(mask))) must match the number of 30-sec epochs ($expected_epochs)."))
        end
    end
    
    # Pre-allocate an array of NamedTuples for type-stable, fast insertion
    results = NamedTuple{(:Epoch, :TimeInSeconds, :Time, :SMax, :AlphaMax, :ThetaMean, :AlphaMean, :SigmaIndex), 
                         Tuple{Int, Float64, Time, Float64, Float64, Float64, Float64, Float64}}[]
    
    for i in eachindex(segs)
        epoch_idx = cld(i, 30)
        
        if mask !== nothing && !mask[epoch_idx]
            continue
        end
        
        amp = AmplitudeSpectrum(segs[i], fs)
        freqs = amp.freq
        spec = amp.spectrum
        
        idx_spindle = findall(f -> 10.5 <= f <= 16.0, freqs)
        idx_alpha_rej = findall(f -> 7.5 <= f <= 10.0, freqs)
        
        S_max = maximum(@view spec[idx_spindle])
        alpha_max = maximum(@view spec[idx_alpha_rej])
        
        if alpha_max > S_max
            continue
        end
        
        idx_theta = findall(f -> 4.0 <= f <= 8.0, freqs)
        idx_alpha = findall(f -> 8.0 <= f <= 10.0, freqs)
        
        theta_mean = mean(@view spec[idx_theta])
        alpha_mean = mean(@view spec[idx_alpha])
        
        val = (2 * S_max) / (alpha_mean + theta_mean)
        
        if val >= threshold
            # Assuming contiguous 1-second segments
            time_start = (i - 1) * 1.0 
            push!(results, (
                Epoch = epoch_idx, 
                TimeInSeconds = time_start, 
                Time = seconds_to_time(time_start),
                SMax = S_max, 
                AlphaMax = alpha_max, 
                ThetaMean = theta_mean, 
                AlphaMean = alpha_mean, 
                SigmaIndex = val
            ))
        end
    end
    
    return DataFrame(results)
end

"""
relative_spindle_power(x::AbstractVector{<:AbstractFloat}, fs::Integer; mask::Union{Nothing, BitVector}=nothing, threshold::Real=0.22)

Computes the Relative Spindle Power (Devuyst et al., 2011) for 1-second segments of EEG data.
Returns a DataFrame containing only the segments where the RSP exceeds the `threshold`.

# Arguments
- `x`: The EEG signal.
- `fs`: Sampling frequency.
- `mask`: Optional BitVector corresponding to 30-sec epochs to include (true) or ignore (false).
- `threshold`: The minimum RSP to be considered a spindle (default: 0.22).

Reference: https://pubmed.ncbi.nlm.nih.gov/22254656/
"""
function relative_spindle_power(ts::TimeSeries; mask::Union{Nothing, BitVector}=nothing, threshold::Real=0.22)
    x = ts.x
    fs = ts.fs
    segs = segment(x, fs; symmetric=true)
    
    if mask !== nothing
        expected_epochs = num_epochs(ts, 30)
        if length(mask) != expected_epochs
            throw(ArgumentError("Mask length ($(length(mask))) must match the number of 30-sec epochs ($expected_epochs)."))
        end
    end
    
    results = NamedTuple{(:Epoch, :TimeInSeconds, :Time, :SpindlePower, :TotalPower, :RSP), 
                         Tuple{Int, Float64, Time, Float64, Float64, Float64}}[]
    
    if isempty(segs)
        return DataFrame(results)
    end

    initial_amp = AmplitudeSpectrum(segs[1], fs, 512)
    freqs = initial_amp.freq
    
    idx_spindle = findall(f -> 11.0 <= f <= 16.0, freqs)
    idx_total = findall(f -> 0.5 <= f <= 40.0, freqs)

    for i in eachindex(segs)
        epoch_idx = cld(i, 30)
        
        if mask !== nothing && !mask[epoch_idx]
            continue
        end
        
        amp = (i == 1) ? initial_amp : AmplitudeSpectrum(segs[i], fs, 512)
        spec = amp.spectrum
        
        spindle_power = sum(@view spec[idx_spindle])
        total_power = sum(@view spec[idx_total])
        
        if total_power > 0
            val = spindle_power / total_power
            if val >= threshold
                time_start = (i - 1) * 1.0
                push!(results, (
                    Epoch = epoch_idx, 
                    TimeInSeconds = time_start, 
                    Time = seconds_to_time(time_start),
                    SpindlePower = spindle_power, 
                    TotalPower = total_power, 
                    RSP = val
                ))
            end
        end
    end
    
    return DataFrame(results)
end
