function compute_tib(stg::Staging; epoch_length=30)
    # Total Time in Bed (TIB) / Total Recording Time
    return length(stg) * epoch_length / 60 # total time in minutes
end

function compute_tst(stg::Staging; epoch_length=30)
    # Total Sleep Time (TST)
    sleep_epochs = count(x -> x in ["1", "2", "3", "4", "5"], stg)
    return sleep_epochs * epoch_length / 60 # sleep time in minutes
end

function compute_sleep_efficiency(stg::Staging; epoch_length=30)
    tst = compute_tst(stg; epoch_length=epoch_length)
    tib = compute_tib(stg; epoch_length=epoch_length)
    
    if tib == 0.0
        return 0.0 # Avoid division by zero
    end
    
    efficiency = (tst / tib) * 100 # sleep efficiency percentage
    return efficiency
end

function compute_sleep_onset(stg::Staging; epoch_length=30)
    i = findfirst(x -> x ∉ ["6", "?"], stg)
    if i === nothing
        throw(ErrorException("No sleep onset found in the staging data."))
    else
        onset_time = (i - 1) * epoch_length / 60  
        return onset_time
    end 
end

function find_waking_epoch(stg::Staging; epoch_length=30)
    N = length(stg)

    for i in N:-1:1
        if stg[i] ∉ ["6", "?"]
            return i + 1
            break
        end
    end
    throw(ErrorException("No waking time found in the staging data."))
end

function compute_wake_after_sleep_onset(stg::Staging; epoch_length=30)
    fell_asleep_epoch = findfirst(x -> x ∉ ["6", "?"], stg)
    # findfirst will return nothing if no sleep onset, but this will be
    # caught by the try/catch in main() via compute_sleep_onset,
    # so we can assume fell_asleep_epoch is valid here.
    
    awoken_epoch = find_waking_epoch(stg; epoch_length=epoch_length)
    waso_epochs = count(x -> x == "6", stg[fell_asleep_epoch:awoken_epoch])
    waso_time = waso_epochs * epoch_length / 60  
    return waso_time
end

function compute_stage_times(stg::Staging; epoch_length=30)
    # Note: We combine N3 and N4 into a single N3/SWS stage
    n1_epochs = count(x -> x == "1", stg)
    n2_epochs = count(x -> x == "2", stg)
    n3_epochs = count(x -> x in ["3", "4"], stg)
    rem_epochs = count(x -> x == "5", stg)

    to_mins(e) = e * epoch_length / 60
    
    return (
        N1_min = to_mins(n1_epochs),
        N2_min = to_mins(n2_epochs),
        N3_min = to_mins(n3_epochs),
        REM_min = to_mins(rem_epochs)
    )
end

function compute_stage_percents(tst::Float64, stage_times::NamedTuple)
    if tst == 0.0
        # Return 0.0 for all percentages if no sleep
        return (N1_pct = 0.0, N2_pct = 0.0, N3_pct = 0.0, REM_pct = 0.0)
    end

    return (
        N1_pct = (stage_times.N1_min / tst) * 100,
        N2_pct = (stage_times.N2_min / tst) * 100,
        N3_pct = (stage_times.N3_min / tst) * 100,
        REM_pct = (stage_times.REM_min / tst) * 100
    )
end

function compute_rem_latency(stg::Staging; epoch_length=30)
    # Time from sleep onset to the first epoch of REM
    # This will throw an error if no onset, which is caught by main()
    sleep_onset_time = compute_sleep_onset(stg; epoch_length=epoch_length)

    i_rem = findfirst(x -> x == "5", stg)
    
    if i_rem === nothing
        return missing # No REM sleep, so REM Latency is missing
    end
    
    rem_start_time = (i_rem - 1) * epoch_length / 60
    
    reml = rem_start_time - sleep_onset_time
    # SOREMP (Sleep Onset REM Period) can result in 0 or negative latency
    # if REM occurs at or before the first defined sleep epoch.
    return reml < 0 ? 0.0 : reml
end

function compute_num_awakenings(stg::Staging)
    # Count transitions to Wake ("6") after sleep onset
    i_onset = findfirst(x -> x ∉ ["6", "?"], stg)
    
    if i_onset === nothing
        return 0 # No sleep onset, so no awakenings after it
    end
    
    count = 0
    for i in (i_onset + 1):length(stg)
        # Check if we *transitioned into* a wake epoch
        if stg[i] == "6" && stg[i-1] != "6"
            count += 1
        end
    end
    return count
end
