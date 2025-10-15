# LyoPronto.jl
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LyoHUB.github.io/LyoPronto.jl/dev)

This package is a Julia complement to [LyoPRONTO](https://github.com/LyoHUB/LyoPronto), an open source Python package.
It has some overlapping functionality with LyoPRONTO, especially simulation of primary drying for conventional lyophilization.
LyoPRONTO (the Python version) also has functionality for generating a design space, estimating time to freeze, and picking optimal drying conditions.
On the other hand, this package has much more advanced utilities for fitting empirical parameters (such as $R_p$ and $K_v$) to experimental data.
This package also provides that fitting functionality for a model applicable to microwave-assisted lyophilization.  

## Installation

From the Julia REPL's Pkg mode (open a REPL and type `]` so that the prompt turns blue), add this package as a Git repo:
```
add https://github.com/LyoHUB/LyoPronto.jl.git
```

## Documentation

The "badge" up above is a link to the documentation, which is [also here](https://lyohub.github.io/LyoPronto.jl/).

## Versioning

In an attempt to adhere to Julia community conventions, this package will use [semantic versioning](semver.org).

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
sol = solve(prob, Rosenbrock23())

modconvtplot(sol)
```

To go beyond one solution to the realm of fitting solutions to experiment, see the docs.