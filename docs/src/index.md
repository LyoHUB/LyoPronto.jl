# Home

## LyoPronto.jl

_A Julia package providing common computations for pharmaceutical lyophilization._

This provides some of the functionality of [LyoPRONTO](https://github.com/LyoHUB/LyoPronto), a Python package.

## Overview

This relatively small package provides a standard literature model for simulating primary drying in pharmaceutical lyophilization, alongside robust utilities for a common parameter estimation workflow. This same infrastructure is provided for a model that adds microwave heating.

Some key advantages this has over the original version of LyoPRONTO are:
- Speed: on my laptop, the regular model can be simulated in about a millisecond. This becomes most relevant when evaluating the model repeatedly in parameter estimation or constructing large design spaces.
- Numerical reliability: in the original LyoPRONTO, bad parameter values (e.g. if input with wrong units or uninformed guesses) easily lead to infinite loops due to the numerical approach used. As a side effect of using a modern library for fast DAE solution, numerical instability errors out instead of hitting an infinite loop. 
- Units: by using `Unitful.jl`, this package enforces dimensional correctness while being compatible with either SI marks or traditional units in lyophilization (like $cm^2\ hr\ Torr / g$ for $R_p$).
- Flexibility: the utilities for fitting $K_v$ and $R_p$ can be used together to fit both at once, not just separately.

## Installation
As a Julia package, this code can be easily installed with the Julia package manager. 

You can add LyoPronto from the Julia General registry (so just like most other packages), using the Julia REPL's Pkg mode (open the REPL and type `]` so the prompt turns blue):
```
add LyoPronto
```
`dev` can be substituted for `add` if you want to make changes to this package yourself, as explained in the [Julia Pkg manual](https://pkgdocs.julialang.org/v1/managing-packages/).

## Dependencies and Reexports

The following are reexported by this package (so that you don't need to import them after importing LyoPronto):
- [OrdinaryDiffEqRosenbrock](https://docs.sciml.ai/DiffEqDocs/stable/); used for solving the DAEs and ODEs inherent here
- [DiffEqCallbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/), used for ending DAE and ODE solves when drying ends
- [Unitful](https://juliaphysics.github.io/Unitful.jl/stable/); specifically, the `u""` macro, `ustrip`, `uconvert`, and `NoUnits`, which is all the API surface needed for regular usage of LyoPronto.

Other noteworthy dependencies:
- [TransformVariables.jl](https://tpapp.github.io/TransformVariables.jl/stable/), which is used to map vector spaces onto realistic parameter values for the inevitable parameter fitting step
- LyoPronto provides [plot recipes](https://docs.juliaplots.org/stable/recipes/) for [Plots.jl](https://docs.juliaplots.org/stable), although it does not depend on Plots in full.  

Plots.jl is used for plotting in this documentation, with the following defaults:
```@example plot_defaults
using Plots # hide
default(:fontfamily, "Computer Modern")
default(:framestyle, :box)
default(:lw, 2)
default(:markersize, 4)
default(:markerstrokewidth, 0.5)
default(:unitformat, :square)
resetfontsizes(); scalefontsizes(1.2)
```



## Authors

Written by Isaac S. Wheeler, a PhD student at Purdue University, advised by Vivek Narsimhan and Alina Alexeenko.
This work was supported in part by funding for NIIMBL project PC4.1-307 .

## Licensing

This package is released with the MIT license.


