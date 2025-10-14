# MicrowaveLyoModeling

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> MicrowaveLyoModeling

To (locally) reproduce this project, do the following:

0. Download this code base. For this project, the data have been included in the 
Git repository, although this is not always done.
1. Open a Julia console in this project directory and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
   This will install all necessary packages for you to be able to run the scripts and
   everything should work out of the box, including correctly finding local paths.

2. Run all the scripts titled `scripts/tuning_[...].jl`. With one or two exceptions, each of these reads in two sets of experimental data, one for conventional and one for microwave-assisted lyophilization, then fits tuning parameters to that experimental data; each script will generate some outputs in `data/exp_pro`, including fit parameters, as well as some figures in `plots`.
3. Run `scripts/comparing_RF_cases.jl`, which generates some comparative plots summarizing all of the experimental cases.


You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate :MicrowaveLyoModeling
```
which auto-activate the project and enable local path handling from DrWatson.

Documentation of `LyoPronto.jl`, which is where the bulk of the math is happening, is available at the following link:

https://lyohub.github.io/LyoPronto.jl/dev/
