# MicrowaveLyoModeling

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.17469773.svg)](https://doi.org/10.5281/zenodo.17469773)

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> MicrowaveLyoModeling
which is the underlying analysis code for a manuscript currently in preparation.

This project is structured as follows:
- The lumped capacitance model is implemented in [LyoPronto.jl](https://lyohub.github.io/LyoPronto.jl/dev/).
- The level set model is implemented in [LevelSetSublimation.jl](https://github.com/Ickaser/LevelSetSublimation.jl), which has LyoPronto.jl as a dependency.
- This "package", in the `src` folder, implements some historical versions of the lumped capacitance model and the model presented in [Srisuma et al., 2023](https://doi.org/10.1016/j.compchemeng.2023.108318).
- This package reexports the following Julia libraries (among others), which are liberally used in the analysis scripts:
  - LyoPronto
  - DrWatson
  - LevelSetSublimation
  - Unitful
  - Plots
  - CSV, TypedTables
  - NonlinearSolve
  - OptimizationOptimJL

To (locally) reproduce this project, do the following:

1. Download this repository, which includes both the data (mostly in `data/exp_raw`) and the analysis code (mostly in the `scripts` folder).
2. Open a Julia console in this project directory and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
   This will install all necessary packages for you to be able to run the scripts and
   everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate :MicrowaveLyoModeling
```
which auto-activate the project and enable local path handling from DrWatson.

If you want to simply rerun all the code, `scripts/rerun_all.jl` will rerun all the figure-generating and analysis scripts.
It runs each script, in order, in its own `let` block (so that each executes in its own local scope). 
This will produce all the graphs in the `plots` folder, as well as putting fit output in `data/exp_pro` and some level set simulations in `data/sims`.

To walk through a particular case, I recommend using [VSCode with the Julia extension](https://code.visualstudio.com/docs/languages/julia) to interactively run a script in order. This will generate some intermediate plots that I do not save for publication, but that are useful for closer inspection.

