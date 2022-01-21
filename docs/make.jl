push!(LOAD_PATH,"../src/")
using Documenter, HammerScalpel

makedocs(modules = [HammerScalpel], sitename = "HammerScalpel.jl")

deploydocs(repo = "github.com/JulianStreyczek/HammerScalpel.jl.git", devbranch = "main")
