

# --------- Constant properties
const Mw = 18.015u"g/mol"
# -------- Approximately constant properties
const ΔHsub = 678u"cal/g"
const θsub = ΔHsub*Mw/Unitful.R

# Thermal conductivity of ice
# From Slack, 1980
# 200K : .032 W/cmK
# 250K : .024 W/cmK
# 273K : .0214 W/cmK
const k_ice = 2.4* u"W/m/K"

const rho_ice = 0.918u"g/cm^3" 
const rho_glass = 2.2u"g/cm^3"
const e_0 = 8.854187e-12u"F/m" # permittivity of free space Coulomb^2/J/m
const σ = 5.670367e-8u"W/m^2/K^4" # Stefan-Boltzmann Constant

# From McCarthy and Fabre, 1989 book chapter
# Thermal conductivity of sucrose, Caster grade powder
const ρ_sucrose = 892u"kg/m^3"
const k_sucrose = 0.139u"W/m/K"

# Water vapor viscosity in dilute limit ---------------
# From Hellmann and Vogel, 2015
# 250K: 8.054 μPa*s
# 260K: 8.383 μPa*s
# 270K: 8.714 μPa*s
const μ = 8.1 * u"μPa*s"

@doc raw"""
    calc_psub(T::F) where F<:Number
    calc_psub(T::Q) where Q<:Quantity

Compute pressure (in Pascals) of sublimation at temperature `T` in Kelvin.

This is essentially an Arrhenius fit, where we compute:
`psub = pref * exp(-ΔHsub*Mw/R/T)`
"""
calc_psub(T::F) where F<:Number = 359.7e10 * exp(-θsub/(T*u"K"))
calc_psub(T::Q) where Q<:Quantity = 359.7e10*u"Pa" * exp(-θsub/uconvert(u"K",T))


"""
Single-purpose module for computing the dielectric loss coefficient of ice, as a function of temperature and frequency.
This module exports a function `ϵpp_f(T, f)` for doing so.
(Set as a separate module just to keep namespaces clean.)
This code is a nearly-direct implementation of the correlation, Eqs. 3-6 presented in:
Takeshi Matsuoka, Shuji Fujita, Shinji Mae; Effect of temperature on dielectric properties of ice in the range 5–39 GHz. J. Appl. Phys. 15 November 1996; 80 (10): 5884–5890. https://doi.org/10.1063/1.363582
"""
module Dielectric
using DataInterpolations
using Unitful
export ϵpp_f

const β = 2.37e4u"K"
const T0 = 15u"K"
const R = 8.3145u"J/mol/K"
const E1 = 55.3u"kJ/mol"
const E2 = 22.6u"kJ/mol"
# const τ0 = T-> (T>223u"K" ? 1.08 : 4.9e7 ) * 5.3e-16u"s" # Fudge factors for high-temp portion of correlation
# const arrhenius = T-> (T > 223u"K" ? exp(E1/R/T) : exp(E2/R/T) )
# const fr = T->1/(2π*τ0(T)*arrhenius(T)) 
# const debye_comp = T-> (β/(T-T0))
# const A = T-> debye_comp(T) * fr(T)
# function ϵpp_f(T, f)
#     uconvert(NoUnits, A(T)/f) + B_interp(T)*ustrip(u"GHz", f)^C_interp(T)
# end

const Tref = [190, 200, 220, 240, 248, 253, 258, 263, 265]*u"K"
const Aref = [0.005, 0.010, 0.031, 0.268, 0.635, 1.059, 1.728, 2.769, 3.326]*1e-4u"GHz"
const Bref = [1.537, 1.747, 2.469, 3.495, 4.006, 4.380, 4.696, 5.277, 5.646]*1e-5
const Cref = [1.175, 1.168, 1.129, 1.088, 1.073, 1.062, 1.056, 1.038, 1.024]
const B_interp = LinearInterpolation(Tref, Bref, extrapolate=true)
const C_interp = LinearInterpolation(Tref, Cref, extrapolate=true)

function ϵpp_f(T, f)
    # The fudge factors of 1.08 and 4.9e7 are my own empirical addition 
    # because the described correlation is off from the values shown in figures
    arrh = 5.3e-16*u"s" * (T > 223u"K" ? 1.08*exp(E1/R/T) : 4.9e7*exp(E2/R/T))
    fr = 1/(2π*arrh)
    A = (β/(T-T0)) * fr
    return uconvert(NoUnits, A/f) + B_interp(T)*ustrip(u"GHz", f)^C_interp(T)
end

end # module Dielectric

# Test that it looks like the figure
# freq = (10 .^ range(7, 11, step=0.1))*u"Hz"
# e1 = ϵpp_f.(200u"K", freq)
# e2 = ϵpp_f.(258u"K", freq)