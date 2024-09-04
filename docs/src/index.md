# Home

## LyoPronto.jl

_A Julia package providing common computations for pharmaceutical lyophilization._

## Overview

This relatively small package puts together some convenience functions useful for simulating primary drying in lyophilization.

## Dependencies and Reexports
It leverages the strengths of the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) ecosystem to do this quickly and efficiently, although it only actually depends on `OrdinaryDiffEqRosenbrock` and `DiffEqCallbacks`, which are both reexported.

Also provided are plot recipes for [Plots.jl](https://docs.juliaplots.org/stable), although this package only depends on `RecipesBase`.

Heavy use is made of [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/), which is reexported.


