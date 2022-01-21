module HammerScalpel

using Distributed, SharedArrays
using NLopt, FiniteDifferences, GridInterpolations
using JLD2, Plots

export solveModel, nopolicy, withpolicy, createFig5, createFig6, createFig7

include("solveModel.jl")
include("nopolicy.jl")
include("withpolicy.jl")
include("withpolicy_parallel.jl")
include("helperfunctions.jl")
include("figures_paper.jl")

end
