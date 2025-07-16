"""
    plot_hypnogram(staging::Vector{String})

Plots a hypnogram from a staging vector. Each element of the vector must be one of:
"1", "2", "3", "4" (NREM stages), "5" (REM), "6" (Wake), or "?" (unstaged).
"""
function plot_hypnogram(staging::Staging)
    # Stage mapping: higher stages plotted lower (hypnogram convention)
    stage_map = Dict("6" => 5,   
                     "5" => 4,   
                     "1" => 3,
                     "2" => 2,
                     "3" => 1,
                     "4" => 0,
                     "?" => NaN) 

    t = 30/(60*60)
    y = map(s -> get(stage_map, string(s), NaN), staging)
    x = collect(0.0:t:t*(length(staging)-1))  # hours if each epoch is 30s

    # Reverse y-axis to show Wake on top
    plot(x, y,
        xlabel = "Time (hours)",
        ylabel = "Stage",
        yticks = ([0, 1, 2, 3, 4, 5], ["N4", "N3", "N2", "N1", "REM", "Wake"]),
        legend = false,
        title = "Hypnogram",
        color=:black,
        lw = 2)
end



