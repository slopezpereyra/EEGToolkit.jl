
"""
    Artifact

A data structure representing a detected artifact in a time series signal.

# Fields

- `loc::Tuple{Int, Int}`: A tuple indicating the start and end sample indices of the artifact within the signal.
- `mean_change::Float64`: A numeric value representing the mean change in signal amplitude during the artifact segment. This quantifies the severity or abruptness of the artifact.
- `test_statistic::Float64`: A test statistic used to determine whether the segment should be considered an artifact.
- `epoch::Integer`: The 30-sec epoch in which the artifact begins.
- `subepoch::Integer`: The 5-sec sub-epoch (within the 30-sec epoch) where the artifact begins.

# Usage

Artifacts are typically detected using automated or manual procedures and stored in a vector of `Artifact` objects. They can be used for visualization, rejection, or correction of contaminated segments in preprocessing pipelines.
"""
struct Artifact
  loc::Tuple{Int, Int}
  mean_change::Float64 
  test_statistic::Float64
  epoch::Integer
  subepoch::Integer
end

const ArtifactData = Union{Vector{Artifact}, Nothing}


"""
function get_epochs_with_artifacts(artifacts::Vector{Artifact})

Given a vector of artifacts, returns a vector of all 30-sec epochs which 
contain an artifact.
"""
function get_epochs_with_artifacts(artifacts::Vector{Artifact})
    epochs = Int[]
    for a in artifacts
        push!(epochs, a.epoch)  
    end
    return sort(unique(epochs))
end



"""
`function plot_artifacts_in_epochs(signal::TimeSeries, artifacts::ArtifactData, 
                                  from::Integer, to::Integer; 
                                  annotate::Bool=false)`

Given a `signal` and a non-empty vector of artifacts, plots the signal and highlights
the existing artifacts from epoch `from` to epoch `to`. If `annotate` is set
to true, artifacts are annotated with their mean change and test statistic.
"""
function plot_artifacts_in_epochs(signal::TimeSeries, artifacts::ArtifactData, 
                                  from::Integer, to::Integer; 
                                  annotate::Bool=false)

    # Vector of artifacts in the desired region
    artifacts = filter(
      x -> x.epoch >= from && x.epoch <= to, 
      artifacts
    )
    if isempty(artifacts)
        @warn "No artifacts found in epoch range"
        return nothing
    end
    
    s, e = epochs_to_indexes(from, to, signal.fs, 30)

    t = (s:e) ./ signal.fs
    sig = signal.x[s:e]

    # Plot base signal without legend and without artifacts
    p = plot(
        t, sig,
        label="",  # no label
        color=:black,
        xlabel="Time (s)",
        ylabel="Amplitude",
        title="Artifacts in Epochs ($from, $to)",
        legend=false,
    )

    # Overlay all artifact regions
    for art in artifacts
        s = art.loc[1]
        e = art.loc[2]

        t_anom = (s:e) ./ signal.fs
        sig_anom = signal.x[s:e]
        plot!(p, t_anom, sig_anom, lw=2, color=:red)

        # Place annotation at the middle of the artifact
        t_mid = mean(t_anom)
        y_mid = maximum(sig_anom)
        
        if annotate
          txt = "Anom \nΔμ=$(round(art.mean_change, digits=1))\nT=$(round(art.test_statistic, digits=1))"
          annotate!(p, t_mid, y_mid, text(txt, :left, 8, :red))
        end
    end

    return p
end


