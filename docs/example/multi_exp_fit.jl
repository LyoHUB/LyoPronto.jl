# # Imports
# LyoPronto is this package. It reexports several other packages, so after
# `using LyoPronto`, you have effectively also done `using Unitful` and a few others.
using LyoPronto

# These packages are used in the test suite,
# but you can use others in their place.

# TypedTables provides a lightweight table structure, not as broadly flexible as a DataFrame but great for our needs.
using TypedTables, CSV
# TransformVariables provides tools for mapping optimization parameters to sensible ranges.
using TransformVariables
# Optimization provides a common interface to a variety of optimization packages, including Optim.
# We import it with OptimizationOptimJL to specify Optim as a backend.
# LineSearches gives a little more granular control over solver algorithms for Optim.
using OptimizationOptimJL
using LineSearches
# Plots is a frontend for several plotting packages, and its companion package StatsPlots has a very nice macro I like. 
using Plots
using StatsPlots: @df
using LaTeXStrings
# Accessors provides the @set and @reset macros for making copies of parameter structs.
using Accessors

# # Get process data
# At this stage, you would want to load in your experimental data, for this example, we will
# generate some synthetic data in order to focus on the fitting approach and assess its
# robustness.
# The scenario we will consider is this: 

# Three different formulations (A, B, and C) are lyophilized in the same lyophilizer,
# with the same vial, fill, and chamber pressure. As a result, $K_v$ should be the same for
# all three experiments, but $R_p$ is likely different for each.


## Vial geometry
## Ran with a 10mL vial, not strictly a 10R but with similar dimensions
ri, ro = get_vial_radii("10R")
Ap = π*ri^2
Av = π*ro^2

## Formulation parameters
csolid = 0.05u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density

## Fill
Vfill = 3u"mL"
hf0 = Vfill / Ap

## Cycle parameters
pch = RampedVariable(100u"mTorr") # constant pressure
T_shelf_0 = -40.0u"°C" # initial shelf temperature
T_shelf_final = -10.0u"°C"  # final shelf temperature
ramp_rate = 0.5 *u"K/minute" # ramp rate
## Ramp for shelf temperature: convert to Kelvin because Celsius doesn't do math very well
Tsh = RampedVariable(uconvert.(u"K", [T_shelf_0, T_shelf_final]), ramp_rate)
## Alternate shelf ramp for third case
TshC = RampedVariable(uconvert.(u"K", [T_shelf_0, T_shelf_final+5u"K"]), 2*ramp_rate)

## For our synthetic cases, we will set values for Kv and Rp
Kv = ConstPhysProp(5.0u"W/m^2/K")
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14.0u"cm*Torr*hr/g"
A2 = 1.0u"1/cm"
RpA = RpFormFit(R0, A1, A2)
RpB = RpFormFit(2R0, 0.5A1, 0.5A2)
RpC = RpFormFit(0.5R0, 2A1, 3A2)

# We will describe each case with a [`ParamObjPikal`](@ref LyoPronto.ParamObjPikal) struct:

poA = ParamObjPikal((
    (RpA, hf0, csolid, ρsolution),
    (Kv, Av, Ap),
    (pch, Tsh)
))
## The @set macro makes a copy of the struct with the specified field changed, 
## so we can use it to make the other two cases.
poB = @set poA.Rp = RpB
poC = @set poA.Rp = RpC
## If we have other parameters to adjust between the cases, we can then use @reset.
@reset poC.Tsh = TshC

# Now, we need to simulate our three synthetic experiments.
sols = [solve(ODEProblem(po), LyoPronto.odealg_chunk2) for po in [poA, poB, poC]]
# Next, load our synthetic data into structs that communicate the fitting problem
t = range(0.0u"hr", stop=maximum([sol.t[end] for sol in sols])*u"hr", step=1u"minute")
fitdats = [PrimaryDryFit(t, (sol.(ustrip.(u"hr", t))[:,1],), t_end=sol.t[end]*u"hr") for sol in sols]

## plot(u"hr", u"°C")
## for (fd, name, color) in zip(fitdats, ["A", "B", "C"], [:blue, :red, :green])
##     plot!(fd)
## end
## plot!(legend=:bottomright)
## savefig("fitdats.svg"); #md #hide
## # ![](fitdats.svg) #md

# # Multi-experiment fitting
# Now, we can set up a multi-experiment fitting problem. The idea is that we will fit the 
# three experiments simultaneously, with a single $K_v$ and three different $R_p$ values. 
# To express this, we will use the `TransformVariables` package to create a transformation 
# that maps a flat vector of numbers to this set of parameters, grouped by experiment.
# LyoPronto provides a function that will apply this transform, then map the parameters to 
# a `ParamObjPikal` for each experiment, and then solve the ODEs and return the objective 
# function value.

## First, a transformation that maps 3 numbers to Rp. LyoPronto provides a convenience function
trans_Rp = Rp_transform_basic(0.75R0, 0.75A1, 0.75A2)

shared_trans = as((separate = as(Vector, trans_Rp, 4),
      shared = as((;Kshf = trans_KRp.Kshf)),
      ))
objnf_pd = OptimizationFunction((x,y)->LyoPronto.objn_pd(x,y,tweight=5e-2), AutoForwardDiff())
