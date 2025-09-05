# # Imports
# LyoPronto is this package. It reexports several other packages, so after
# `using LyoPronto`, you have effectively also done `using Unitful` and a few others.
using LyoPronto

# These are other packages that I use in the test suite,
# but you can use others in their place.
# TypedTables provides a lightweight table structure, not as broadly flexible as a DataFrame but great for our needs
using TypedTables, CSV
# Plots is a frontend for several plotting packages, and its companion package StatsPlots has a very nice macro I like. 
using Plots
using StatsPlots: @df
using LaTeXStrings
# For dealing with parameter structs and making copies, Accessors provides the @set and @reset macros
using Accessors

# # Example Process Data

## Data start at 8th row of CSV file.
## This needs to point to the right file, which for documentation is kinda wonky
procdata_raw = CSV.read(joinpath(@__DIR__, "..", "..", "example", "2024-06-04-10_MFD_AH.csv"), Table, header=8)
t = uconvert.(u"hr", procdata_raw.CycleTime .- procdata_raw.CycleTime[1])
## At midnight, timestamps revert to zero, so catch that case
for i in eachindex(t)[begin+1:end]
    if t[i] < t[i-1]
        t[i:end] .+= 24u"hr"
    end
end
## Some of the dispatches don't like if time is not a float
t = float.(t)

## Rename the columns we will use, and add units
procdata = map(procdata_raw) do row
    ## In the anonymous `do` function, `row` is a row of the table.
    ## Return a new row as a NamedTuple
    (pirani = row.VacPirani * u"mTorr",
     cm = row.VacCPM * u"mTorr",
     T1 = row.TP1 * u"°C",
     T2 = row.TP2 * u"°C",
     T3 = row.TP4 * u"°C", # Quirk of this experimental run: TP3 slot was empty
     Tsh = row.ShelfSetPT * u"°C",
     phase = row.Phase # identify whether freezing, primary drying, or secondary
    )
end
procdata = Table(procdata, (;t)) # Append time to table

## Count time from the beginning of experiment
pd_data = filter(row->row.phase == 4, procdata)
tstart_pd = pd_data.t[1]
pd_data.t .-= pd_data.t[1]

t_end = identify_pd_end(pd_data.t, pd_data.pirani, :onoff)

# # Example models
# ## Conventional lyophilization

## Vial geometry
## Ran with a 10mL vial, not strictly a 10R but with similar dimensions
Ap, Av = π.*get_vial_radii("10R") .^ 2

## Experimental conditions
T_shelf_0 = -40.0u"°C" |> u"K" # initial shelf temperature, in Kelvin for math reasons
T_shelf_final = -10.0u"°C" |> u"K" # final shelf temperature
ramp_rate = 0.5 *u"K/minute" # ramp rate
## Set points, followed, by ramp rate, followed by hold times if there are multiple ramps
Tsh = RampedVariable([T_shelf_0, T_shelf_final], ramp_rate)
## Single set point with no ramps
pch = RampedVariable(100u"mTorr")

## Formulation parameters
csolid = 0.05u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density
## Previously fitted values for Rp
R0 = 0.93u"cm^2*Torr*hr/g"
A1 = 21.1u"cm*Torr*hr/g"
A2 = 1.2u"1/cm"
Rp = RpFormFit(R0, A1, A2)
## Fit value for heat transfer coeff
Kshf = ConstPhysProp(13.9u"W/m^2/K")

## Fill
Vfill = 3u"mL"
hf0 = Vfill / Ap

po = ParamObjPikal([
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh)
]);

prob = ODEProblem(po)
sol_conv = solve(prob, Rodas3());

# # System setpoints and conditions with `RampedVariable`s

# Because it is very common to have a gradual ramp in temperature at the start of drying 
# (and in general any time the set point changes), LyoPronto provides a tool for concisely 
# describing the set point over time.

# For a constant set point (with no ramps), provide a single value:
Tsh1 = RampedVariable(-15u"°C")
# For a ramp from freezing temperature followed by a single hold (as is common in primary 
# drying), provide the initial temperature, target temperature, and ramp rate:
## Convert from Celsius to Kelvin so that algebra can happen on backend
Tsh2 = RampedVariable([-40.0, -10]u"°C" .|> u"K", 1.0u"K/minute")

# For multiple ramps, provide set points, then ramp rates, then hold times in between ramps.
# In general, supply one less ramp than set points, and one less hold time than ramps.
Tsh3 = RampedVariable([-40, -20, 0]u"°C" .|> u"K", [2//3, 0.5]u"K/minute", [1u"hr"])

# A plot recipe is provided for all of these `RampedVariable`s.
plot(xunit=u"hr", xlimit=(-0.1u"hr", 2.5u"hr"))
plot!(Tsh1, lw=3, label="No ramp")
plot!(Tsh2, lw=3, label="1 ramp")
plot!(Tsh3, lw=3, label="2 ramps")

# # Estimating Rp over time

# The standard (Pikal) model for primary drying in lyophilization consists of, fundamentally,
# one ODE (change in frozen layer height) with a nonlinear algebraic constraint (pseudosteady 
# heat and mass transfer, coupled at sublimation front).
# Since we have a system with one ODE and one nonlinear algebraic constraint, there is only 
# one degree of freedom at a given time point. So with temperature measurements over time, 
# we can compute corresponding mass flow over time or mass transfer resistance over length.

# In LyoPronto, this functionality is implemented with the `[calc_hRp_T](@ref)` function,
# (think "compute $Rp(h_d)$ from $T_f(t)$").

# To do so, we need to know about the experimental conditions; for that purpose, we pass a 
# `ParamObjPikal` containing that information.
# To deal with the actual temperature series, we use a `PrimaryDryFit` object, which allows
# us to encode the way that, at some point, each temperature series deviates from the 
# regular pseudosteady behavior governed by this model.
# (Strictly speaking, there are a variety of phenomena involved, but for here it is enough 
# to say that at some point in time each temperature series experiences a sharp rise that is
# not described by the model.)

## The PrimaryDryFit:
fitdat_all = @df pd_data PrimaryDryFit(:t, (:T1[:t .< 15u"hr"],
                                    :T2[:t .< 13u"hr"],
                                    :T3[:t .< 16u"hr"]),)
plot(fitdat_all)

# Note that in this plot, T1 rises after 13 hours--I have deliberately included that
# to show what this will do in $R_p(h_d)$ space.
# Now, with the `ParamObjPikal` and `PrimaryDryFit` defined, we can calculate $R_p(h_d)$:

## Compute just for the first temperature series
calc_hRp_T(po, fitdat_all, i=1)
## Compute for all temperature series
hRps = [calc_hRp_T(po, fitdat_all; i) for i in 1:3]

# This returned a vector tuples, with a vector each for `h_d(t)` and `R_p(t)` for each set
# of temperature measurements. To plot this against a temperature fit as in the other 
# tutorial here, we can do the following:

pl = plot(xlabel="h_d", ylabel="R_p", xunit=u"cm", yunit=u"cm^2*Torr*hr/g")
for (i, hRp) in enumerate(hRps)
    plot!(hRp[1], hRp[2], label="T$i")
end
## For comparison, plot Rp as computed from fit in the other example
l = range(0u"cm", hf0, length=100)
plot!(l, Rp.(l), label="Direct fit to \$T(t)\$")
