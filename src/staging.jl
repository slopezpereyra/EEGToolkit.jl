
"""
Structure for staging data. Essentially, any vector `x` of strings satisfying 

        Set(unique(x)) ⊈ Set(["1", "2", "3", "4", "5", "6", "?"])

or 

        v ∈ x ⟹  v ∈ {1, 2, 3, 4, 5, 6, ?}

where the set is a set of strings. By convention, in a staging vector:

- "n" is the nth stage of sleep for 1 ≤ n ≤ 4 
- 5 is REM 
- 6 is wakefulness.
- The ith element in the vector corresponds to the stage of the ith 30-sec 
epoch (in some EEG).

Constructor: Staging(stages::Vector{String})
"""
struct Staging <:AbstractVector{String}
  stages::Vector{String}

  function Staging(stages::Vector{String})

    if Set(unique(stages)) ⊈ Set(["1", "2", "3", "4", "5", "6", "?"])
      e = "Invalid stage score: Sleep stages should be 1, 2, 3, 4, 5, 6, ?, 
        where 5 marks REM, 6 wakefulness, and ? unknown/unscored."
      throw( ArgumentError(e) )
    end
    new(stages)
  end
end

Base.size(s::Staging) = size(s.stages)
Base.getindex(s::Staging, i::Int) = s.stages[i]
Base.setindex!(s::Staging, v::String, i::Int) = (s.stages[i] = v)
Base.length(s::Staging) = length(s.stages)
Base.iterate(s::Staging, state...) = iterate(s.stages, state...)

"""
`function transition_matrix(s::Staging)`

Compute the transition matrix of a staging vector.
"""
function transition_matrix(s::Staging;
                           from_labels=["6", "5", "1", "2", "3"],
                           stage_names=["AWA", "REM", "N1", "N2", "N3"])

    idx = Dict(stage => i for (i, stage) in enumerate(from_labels))
    mat = zeros(Int, length(from_labels), length(from_labels))

    for i in 1:(length(s)-1)
        from = s[i]
        to = s[i+1]
        if haskey(idx, from) && haskey(idx, to)
            mat[idx[from], idx[to]] += 1
        end
    end

    return mat, stage_names
end

"""
`function plot_transition_matrix(s::Staging)`

Plots the transition matrix of a staging vector.
"""
function plot_transition_matrix(s::Staging)
    mat, labels = transition_matrix(s)
    n = length(labels)

    heatmap(
        labels, labels, mat;
        xlabel="To",
        ylabel="From",
        title="Sleep Stage Transition Matrix",
        c=:blues,
        size=(600, 500),
        right_margin=5Plots.mm,
        colorbar_title="Count",
        aspect_ratio=1,
        yflip=true
    )
    print(mat)
    annotate!(
        vec(tuple.((1:n)'.-0.5, (1:n).-0.5, string.(mat)))
    )
        
end
