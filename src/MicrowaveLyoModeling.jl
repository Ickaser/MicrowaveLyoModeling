module MicrowaveLyoModeling

using Reexport
@reexport using LyoPronto
@reexport using LevelSetSublimation
@reexport using DrWatson
@reexport using Unitful
@reexport using Plots
@reexport using StatsPlots: @df
@reexport using LaTeXStrings
@reexport using CSV
@reexport using TypedTables
@reexport using NonlinearSolve
@reexport using OptimizationOptimJL
@reexport using SavitzkyGolay
@reexport using LaTeXStrings
@reexport using DataInterpolations: LinearInterpolation, ExtrapolationType, ConstantInterpolation
@reexport using Accessors
@reexport using PrettyTables
using PrecompileTools
using SpecialFunctions: besselj0, besselj1
using Roots
using Dates

using LineSearches
optalg = LBFGS(linesearch=LineSearches.BackTracking())
export optalg
const objf_KBB = OptimizationFunction(LyoPronto.obj_pd, AutoForwardDiff(chunksize=3))
const objf_Rp = OptimizationFunction(LyoPronto.obj_pd, AutoForwardDiff(chunksize=3))
const objf_KRp = OptimizationFunction(LyoPronto.obj_pd, AutoForwardDiff(chunksize=4))
export objf_KBB, objf_Rp, objf_KRp

include(srcdir("plotrecipes.jl"))
include(srcdir("braatz.jl"))
include(srcdir("LC_alternate.jl"))
include(srcdir("exp_data_str.jl"))
# include(srcdir("precompilation.jl"))

end