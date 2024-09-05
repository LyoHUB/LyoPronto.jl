# Home

## LyoPronto.jl

_A Julia package providing common computations for pharmaceutical lyophilization._


## Overview

This relatively small package puts together some convenience functions useful for simulating primary drying in lyophilization.

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


