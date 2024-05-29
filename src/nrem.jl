include("eeg.jl")

function insert_spaces_on_change(str)
    # Initialize an empty string to build the modified version
    modified_str = ""

    # Iterate through the characters of the input string
    for (index, char) in enumerate(str)
        # Append the current character to the modified string
        modified_str *= char

        # Check if we're at the last character or the next character is different
        if index < length(str) && str[index] != str[index+1]
            # Append a space to separate different character sequences
            modified_str *= " "
        end
    end

    return modified_str
end

function stage_to_word(stages)

    # Given a stages s₁, ..., sₙ, produces a sentence s.t. each word 
    # in a sentence corresponds to segments of the staging with repeating 
    # stages.
    stages = replace(stages, "1" => "I", "?" => "I", "2" => "X", "3" => "X", "4" => "X")
    α = insert_spaces_on_change(join(stages, ""))

    # A function to transform the words of α, effectively mapping α 
    # onto a language more suited for regexes.
    function transform(w)
        if occursin('I', w) || occursin("X", w)
            return w
        end
        # If previous condition was not met, we have either `5` or `6`.
        if length(w) >= 10
            return (w)
        end
        repeat("I", length(w))
    end

    trans = transform.(split(α))
    return join(trans, "")
end


function nrem(staging::Vector, n::Integer=30, m::Integer=10)

    if Set(unique(staging)) ⊈ Set(["1", "2", "3", "4", "5", "6", "?"])
        @error "Invalid stage score: Sleep stages should be 1, 2, 3, 4, 5, 6, ?, 
        where 5 marks REM, 6 wakefulness, and ? unknown/unscored."
    end

    α = stage_to_word(staging)
    φ = Regex("(?:[X](I*)){$n,}(?:5|6{$m,}|\$)")
    ζ = Regex("(?:[X](I*)){$n,}(?:[5](I*){$m,}|6{$m,}|\$)")
 
    # ------------------------------------
    # Finding first and last NREM periods:
    x = eachmatch(φ, α)
    locs = []
    for match in x
        s = match.offset
        e = match.offset + length(match.match) - 1
        push!(locs, [s, e])
    end

    # Keep only first and last matches.
    locs = [locs[1], locs[length(locs)]]
    s = (locs[1][2])
    e = locs[length(locs)][1] - 1
    
    # ------------------------------------
    # Finding middle NREM periods.

    # Split α so as to search for ζ in between the first and the last matches
    middle = α[s:e]
    matches = eachmatch(ζ, middle)

    for match in matches
        begins = s + match.offset - 1
        loc = [begins, begins + length(match.match)]
        insert!(locs, length(locs), loc)
    end

    # Collect all matches into vectors
    locs = [collect(x[1]:x[2]) for x in locs]

    # Remove ending REM and Wakes as well as in-between REM, Wake, and Stage 1 epochs.
    for nrem in locs
        indx = findall(x -> staging[x] != "2" && staging[x] != "3" && staging[x] != "4", nrem)
        deleteat!(nrem, indx)
    end

    return locs
end

function nrem(eeg::EEG, n::Integer=30, m::Integer=10)
    nrem(eeg.staging, n, m)
end
