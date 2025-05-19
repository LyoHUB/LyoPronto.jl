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

From the Julia REPL's Pkg mode (open a REPL and type `]` so that the prompt turns blue), add this package as a Git repo:
```
add https://github.com/LyoHUB/LyoPronto.jl.git
```
`dev` can be substituted for `add` if you want to make changes to this package yourself, as explained in the [Julia Pkg manual](https://pkgdocs.julialang.org/v1/managing-packages/).

## Dependencies and Reexports
This package leverages the strengths of the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) ecosystem to solve equations quickly and efficiently, although it only directly depends on `OrdinaryDiffEqRosenbrock` and `DiffEqCallbacks`, which are both reexported.

Also provided are plot recipes for [Plots.jl](https://docs.juliaplots.org/stable), although this package only depends on `RecipesBase`.

Heavy use is made of [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/), which is reexported.


## Authors

Written by Isaac S. Wheeler, a PhD student at Purdue University.
This work was supported in part by funding for NIIMBL project PC4.1-307 .

## License

None yet. My intentions are to use the MIT license once this has been published in a scientific journal.


