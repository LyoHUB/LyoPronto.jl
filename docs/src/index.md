# Home

## LyoPronto.jl

_A Julia package providing common computations for pharmaceutical lyophilization._

This package is a Julia complement to [LyoPRONTO](https://github.com/LyoHUB/LyoPronto), a Python package [shivkumarLyoPRONTOOpenSourceLyophilization2019](@cite) with a [web interface](https://lyopronto.geddes.rcac.purdue.edu). It is not a one-to-one translation, but rather a reimplementation of the same underlying model, with a much-improved interface for parameter estimation and an interface for extending that infrastructure to new models.

In the newer [web interface](https://lyopronto2.geddes.rcac.purdue.edu), some of the functionality is provided by calling Python and some is provided by this package.

## Overview

Some key advantages this has over the original (Python) version of LyoPRONTO are:
- Speed: on my laptop, the regular model can be simulated in about a millisecond. This becomes most relevant when evaluating the model repeatedly in parameter estimation or constructing large design spaces (both of which take less than a second for a well-posed problem).
- Numerical reliability: This version uses `OrdinaryDiffEq.jl` for solving the ODEs and DAEs, which is a modern and robust library for fast numerical solution. This provides a lot of bells and whistles which we actively use, on top of being performant. 
- Units: by using `Unitful.jl`, this package enforces dimensional correctness while being compatible with either SI marks or traditional units in lyophilization (like $cm^2\ hr\ Torr / g$ for $R_p$).
- Flexibility: the utilities for fitting parameters like $K_v$ and $R_p$ can be used together to fit both at once, not just separately, and temperature data can be used in conjunction with drying time data to constrain rigorous least-squares fits.
- Extensibility: the package provides a framework for implementing new physical models for processes similar to freeze drying, and defining new experimental data types for parameter fitting.

As a consequence (and motivating example) of its extensibility, LyoPronto.jl also implements a model for radio frequency-assisted lyophilization, which is not available in the original LyoPRONTO.

## Installation
As a Julia package, this code can be easily installed with the Julia package manager. 

You can add LyoPronto from the Julia General registry (so just like most other packages), using the Julia REPL's Pkg mode (open the REPL and type `]` so the prompt turns blue):
```
add LyoPronto
```
`dev` can be substituted for `add` if you want to make changes to this package yourself, as explained in the [Julia Pkg manual](https://pkgdocs.julialang.org/v1/managing-packages/).

## Dependencies and Reexports

Among the dependencies of LyoPronto are a few packages which provide functionality without which LyoPronto would be unusable, so those functions are exported by LyoPronto as well (so that `using LyoPronto` makes these functions available). This includes the following:
- From [OrdinaryDiffEqRosenbrock](https://docs.sciml.ai/DiffEqDocs/stable/) and [OrdinaryDiffEqNonlinearSolve](https://docs.sciml.ai/DiffEqDocs/stable/), `ODEProblem`, `solve` used for solving the DAEs and ODEs inherent here. The recommended algorithm for LyoPronto's systems is now `Rodas5P(AutoForwardDiff(chunksize=2))`, which is made public as `LyoPronto.odealg_chunk2`.
- [Unitful](https://juliaphysics.github.io/Unitful.jl/stable/); specifically, the `u""` macro, `ustrip`, `uconvert`, and `NoUnits`, which is all the API surface needed for regular usage of LyoPronto.

Other noteworthy dependencies:
- [TransformVariables.jl](https://tpapp.github.io/TransformVariables.jl/stable/), which is used to map vector spaces onto realistic parameter values for the inevitable parameter fitting step
- LyoPronto provides [plot recipes](https://docs.juliaplots.org/stable/recipes/) for [Plots.jl](https://docs.juliaplots.org/stable), although it does not depend on Plots in full.  

Plots.jl is used for plotting in this documentation.

## Authors

Written by Isaac S. Wheeler, a PhD student at Purdue University, advised by Vivek Narsimhan and Alina Alexeenko.
This work was supported in part by funding for NIIMBL project PC4.1-307 .

## Licensing

This package is released with the MIT license.

## Cited References

The following works are cited in this documentation
```@bibliography
```
