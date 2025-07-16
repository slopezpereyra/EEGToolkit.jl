
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
