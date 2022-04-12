module Delta2
using DataFrames
using CSV
using RecipesBase
using Measures
import REPL
using REPL.TerminalMenus
import StatsBase.geomean
import Statistics.mean
include("ddct.jl")
include("plotting.jl")
include("wizards.jl")
export QPCRDataset,DeltaCT, DDCT, targets, samples, writeresult, readpcr
end # module
