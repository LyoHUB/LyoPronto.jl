# LyoPronto.jl
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LyoHUB.github.io/LyoPronto.jl/dev)
[![](https://zenodo.org/badge/DOI/10.5281/zenodo.17373491.svg)](http://doi.org/10.5281/zenodo.17373491)

_A Julia package providing common computations for pharmaceutical lyophilization._

This package is a Julia complement to [LyoPRONTO](https://github.com/LyoHUB/LyoPronto), a Python package with a [web interface](https://lyopronto.geddes.rcac.purdue.edu). It is not a one-to-one translation, but began as a reimplementation of the same underlying mathematical model, with a much-improved interface for parameter estimation and an interface for extending that infrastructure to new models.

In the newer [web interface](https://lyopronto2.geddes.rcac.purdue.edu), some of the functionality is provided by calling Python and some is provided by this package.

## Overview

It has some overlapping functionality with LyoPRONTO, especially simulation of primary drying for conventional lyophilization.
LyoPRONTO (the Python version) also has functionality for generating a design space, estimating time to freeze, and picking optimal drying conditions.

Some key advantages this has over the original (Python) version of LyoPRONTO are:
- Speed: on my laptop, the regular model can be simulated in about a millisecond. This becomes most relevant when evaluating the model repeatedly in parameter estimation or constructing large design spaces (both of which take less than a second for a well-posed problem).
- Numerical reliability: This version uses `OrdinaryDiffEq.jl` for solving the ODEs and DAEs, which is a modern and robust library for fast numerical solution. This provides a lot of bells and whistles which we actively use, on top of being performant. 
- Units: by using `Unitful.jl`, this package enforces dimensional correctness while being compatible with either SI marks or traditional units in lyophilization (like $cm^2\ hr\ Torr / g$ for $R_p$).
- Flexibility: the utilities for fitting parameters like $K_v$ and $R_p$ can be used together to fit both at once, not just separately, and temperature data can be used in conjunction with drying time data to constrain rigorous least-squares fits.
- Extensibility: the package provides a framework for implementing new physical models for processes similar to freeze drying, and defining new experimental data types for parameter fitting.

As a consequence (and motivating example) of its extensibility, LyoPronto.jl also implements a model for radio frequency-assisted lyophilization, which is not available in the original LyoPRONTO.

## Installation

As a Julia package, this code can be easily installed with the Julia package manager. 

From the Julia REPL's Pkg mode (open a REPL and type `]` so that the prompt turns blue), add this package from the General registry with:
```
add LyoPronto
```


## Documentation

The "badge" up above is a link to the documentation, which is [also here](https://lyohub.github.io/LyoPronto.jl/).

## Versioning

In accordance with the Julia community's conventions, this package uses [semantic versioning](semver.org).

## Authors

Written by Isaac S. Wheeler, a PhD student at Purdue University, advised by Prof. Vivek Narsimhan and Prof. Alina Alexeenko. 
This work was supported in part by funding for NIIMBL project PC4.1-307 .

## License

MIT License; see `LICENSE` file.

# Example usage

```julia
using LyoPronto

# Vial information
Ap, Av = @. π*get_vial_radii("6R")^2  # cross-sectional area inside the vial
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)

# Formulation parameters
csolid = 0.06u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g" # Guess
A1 = 14.0u"cm*Torr*hr/g" # Guess
A2 = 1.0u"1/cm" # Guess
Rp = RpFormFit(R0, A1, A2)

# Cycle parameters
Vfill = 3u"mL" # ml
pch = RampedVariable(70u"mTorr")
Tsh = RampedVariable([-15u"°C", 10u"°C"].|>u"K", 0.5u"K/minute")
hf0 = Vfill / Ap

# Put information together
po = ParamObjPikal((
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh)
))

prob = ODEProblem(po)
sol = solve(prob, LyoPronto.odealg_chunk2)

modconvtplot(sol)
```

To go beyond one solution to the realm of fitting solutions to experiment, see the docs.