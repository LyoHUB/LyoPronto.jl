# Example: Tuning Mass Transfer Resistance

The source code for this example can be found in `LyoPronto.jl/scripts/tuning_sucrose_Rp.jl`. 

## Setup

My current recommended approach for making use of this package is to have installed `LyoPronto` as a dependency of your project, according to the [install instructions](@ref "Installation")

For management of your project, [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) is an effective tool which I like. I will demonstrate some usage of that tool here.

## Packages Used

The following packages are used in this example:

First, load packages:
```julia
using DrWatson
using LyoPronto
using CSV
using Plots
using LaTeXStrings
using SavitzkyGolay
using Optim
```

## Loading Experimental Data

First, load a `.csv` of experimental data with the `CSV.File` method, and use `DrWatson`'s `datadir` function to help keep things clean.

This file has 7 rows of information at the top, so the data headers are at row 8; it is located in the project directory under `[project]/data/exp_raw/"2024-06-21-16_MFD_AH.csv`.

```julia
procdat = CSV.File(datadir("exp_raw", "2024-06-21-16_MFD_AH.csv"), header=8)
```

For the Millrock MicroFD lyophilizer that data was gathered on, the `.csv` file has a column `Phase` indicating what step of the process is going on, so to get our primary drying data out we can identify it by rows where the `Phase == 4`:

```julia
is_PD = procdat["Phase"] .== 4
```

To get the data series, we index to those rows for each variable of interest, and mark their units from `Unitful.jl` by multiplying by `u"[unit]"`.
The time points are once a minute, so here we construct a set of time data from scratch rather than mess with reading in the `DateTime` values from the lyophilizer (which may flip to zero at midnight and have other annoying behavior).

```julia
tproc = range(0, length=sum(is_PD), step=1/60)*u"hr"

p_pir = procdat["VacPirani"][is_PD]*u"mTorr"
Tsh_d = procdat["ShelfSetPT"][is_PD]*u"°C"
T1 = procdat["TP1"][is_PD]*u"°C"
T2 = procdat["TP2"][is_PD]*u"°C"
T3 = procdat["TP4"][is_PD]*u"°C"
```

To check that everything looks right, plot the temperatures, taking advantage of a recipe from this package, as well as the `L"[latex]"` macro from `LaTeXStrings`. 
(This requires `Plots.jl`.)

```julia
tplotexperimental(tproc, T1, T2, T3, labels=[L"T_{p1}" L"T_{p2}" L"T_{p3}"])
plot!(tproc, Tsh_d, c=:black)
```

Based on an examination of the temperature data, we want to go only up to the "temperature rise" commonly observed in lyophilization near (but not at) the end of drying. This happens at 7.5 hours for `T1` and `T3` and at about 12.5 hours for `T2`.

```julia
T1trm = T1[tproc .< 7.5u"hr"]
T2trm = T2[tproc .< 12.5u"hr"]
T3trm = T3[tproc .< 7.5u"hr"]
```


Look at the Pirani pressure and ascertain the end of drying by using a Savitzky-Golay filter to identify a maximum in the second derivative.
Note that because `SavitzkyGolay` doesn't play nice with Unitful, we strip out units and add them back in.
Another plot recipe plots the end of drying as a vertical line on that pressure graph.

```julia
p_pir_sm = savitzky_golay(ustrip.(u"mTorr", p_pir), 91, 3, deriv=0).y *u"mTorr" # 91 a window width; 3 the polynomial order
p_pir_der2 = savitzky_golay(ustrip.(u"mTorr", p_pir), 91, 3, deriv=2).y
t_end = tproc[argmax(p_pir_der2[100:end])+99]
plot(tproc, p_pir, label="data")
plot!(tproc, p_pir_sm, label="smoothed data")
tendplot!(t_end, label=L"t_{end}")
```

## Set up model parameters

Below, we make liberal use of Unitful units. Also note that [`RampedVariable`](@ref) and [`RpFormFit`](@ref LyoPronto.RpFormFit) are used to simplify some common things we use.
```julia
# Vial geometry
vialsize = "6R"
Ap, Av = π.*get_vial_radii(vialsize).^2

# Formulation parameters; Rp here is a placeholder guess
c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14u"cm*Torr*hr/g"
A2 = 1u"1/cm"
Rp = RpFormFit(R0, A1, A2)

# Fill
Vfill = 3u"mL" # ml
hf0 = Vfill / Ap

# Cycle parameters
pch = RampedVariable(70u"mTorr") # constant pressure
T_shelf_0 = -40.0u"°C" # initial shelf temperature
T_shelf_final = -10.0u"°C"  # final shelf temperature
ramp_rate = 0.5 *u"K/minute" # ramp rate
# Ramp for shelf temperature: convert to Kelvin because Celsius doesn't do math very well
Tsh = RampedVariable(uconvert.(u"K", [T_shelf_0, T_shelf_final]), ramp_rate)

# Guess for heat transfer
KC = 6.556e-5u"cal/s/K/cm^2"
KP = 2.41e-3u"cal/s/K/cm^2/Torr"
KD = 2.62u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
# Alternative guess, without pressure dependence:
# Kshf = p-> 5.0u"W/m^2/K"

params_bunch = [
    (Rp, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh)
]

```

## Run a sanity-check simulation

```julia
# Time span: used to set initial time and to give an upper bound on time, in case parameters are bad
tspan = (0.0, 100.0) # hours
# Initial condition
u0 = [ustrip(u"cm", hf0), ustrip(u"K", Tsh(0u"minute"))]

# Set up as an ODE problem
prob = ODEProblem(lyo_1d_dae_f, u0, tspan, tuple(params_bunch...))
# Solve with the Rodas4P() algorithm, use a callback to terminate at end of drying
sol = solve(prob, Rodas4P(), callback=LyoPronto.end_drying_callback)

```
And then plot these results with a recipe, to make sure that the chosen $K_{sh-f}$ and $R_p$ are sane starting points:

```julia
plot(tproc, Tsh_d, c=:black, label=L"T_{sh}")
tplotmodelconv!(sol, label=L"$T_p$, model")
```

## Minimize least square difference to compare to experimental data

First, set up a function that takes $K_{sh-f}$ and $R_p$ and returns a solution object (see [`gen_sol_conv_dim`](@ref LyoPronto.gen_sol_conv_dim) ).

```julia
otherparams = (hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)
gen_sol_conv = KRp -> gen_sol_conv_dim(KRp, otherparams, u0, tspan)
```

Next, set up an objective function we will minimize (see [`obj_tT_conv`](@ref LyoPronto.obj_tT_conv) for details).
We will compare to `T1`.

```julia
obj_KRp1 = KRp -> obj_tT_conv(KRp, gen_sol_conv, (tproc[tproc.<7.5u"hr"], T1trm), t_end=t_end) 
```

Minimize this objective function, using `Optim` with the `NelderMead()` algorithm:
```julia
opt_KRp1 = optimize(obj_KRp1, [12, 0.1, 5, 0.1], NelderMead())
```

Get out the found values of our tuning parameters, and generate the corresponding solution profile:

```julia
KRp_1 = Optim.minimizer(opt_Krp1)
prof1 = gen_sol_conv(Krp_1)[1]
```

Plot and compare experiment to model:
```julia
tplotexperimental(tproc, T1, T2, T3)
tplotmodelconv!(prof1)
tendplot!(t_end)
```


