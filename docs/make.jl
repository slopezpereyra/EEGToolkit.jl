# Inside make.jl
push!(LOAD_PATH,"../src/")
using EEGToolkit
using Documenter

makedocs(
         sitename = "EEGToolkit.jl",
         modules  = [EEGToolkit],
         pages=[
                "Home" => "index.md"
               ],
        clean = true)
deploydocs(;
    repo="github.com/slopezpereyra/EEGToolkit.jl.git",
)
