# LyoPRONTO.jl
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LyoHUB.github.io/LyoPronto.jl/dev)

This package is a Julia reimplementation of some core functionality of [LyoPRONTO](https://github.com/LyoHUB/LyoPronto), an open source Python package.

## Installation

From the Julia REPL's Pkg mode (open a REPL and type `]` so that the prompt turns blue), add this package as a Git repo:
```
add https://github.com/LyoHUB/LyoPronto.jl.git
```

## Documentation

The "badge" up above is a link to the documentation.

## Versioning

In an attempt to adhere to Julia community conventions, this package will use [semantic versioning](semver.org).

## Authors

Written by Isaac S. Wheeler, a PhD student at Purdue University.
This work was supported in part by funding for NIIMBL project PC4.1-307 .

## License

None yet. My intentions are to use the MIT license once this has been published in a scientific journal.

# Example usage

```julia
using Optimization, OptimizationOptimJL
using LineSearches
using NonlinearSolve
optalg = LBFGS(linesearch=LineSearches.BackTracking())

# Vial information
Ap, Av = @. π*get_vial_radii("6R")^2  # cross-sectional area inside the vial
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)

# Formulation parameters
csolid = 0.06u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14.0u"cm*Torr*hr/g"
A2 = 1.0u"1/cm"
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



```