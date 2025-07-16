"""
    plot_hypnogram(staging::Vector{String})

Plots a hypnogram from a staging vector. Each element of the vector must be one of:
"1", "2", "3", "4" (NREM stages), "5" (REM), "6" (Wake), or "?" (unstaged).
"""
function plot_hypnogram(staging::Staging)
    # Stage mapping: higher stages plotted lower (hypnogram convention)
    stage_map = Dict("6" => 0,   # Wake
                     "5" => 1,   # REM
                     "4" => 2,   # N3/N4
                     "3" => 3,   # N3
                     "2" => 4,   # N2
                     "1" => 5,   # N1
                     "?" => NaN) # Unstaged

    t = 30/(60*60)
    y = map(s -> get(stage_map, string(s), NaN), staging)
    x = collect(0.0:t:t*(length(staging)-1))  # hours if each epoch is 30s

    # Reverse y-axis to show Wake on top
    plot(x, y,
        xlabel = "Time (hours)",
        ylabel = "Stage",
        yticks = ([0, 1, 2, 3, 4, 5], ["Wake", "REM", "N4", "N3", "N2", "N1"]),
        legend = false,
        title = "Hypnogram",
        color=:black,
        lw = 2)
end

