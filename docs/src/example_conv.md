# Example: Tuning Mass Transfer Resistance

The source code for this example can be found in `LyoPronto.jl/scripts/tuning_sucrose_Rp.jl`. 

## Setup

My current recommended approach for making use of this package is to have installed `LyoPronto` as a dependency of your project, according to the [install instructions](@ref "Installation")

For management of your project, [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) is an effective tool which I like. I will demonstrate some usage of that tool here.

## Packages Used

First, we will load packages.
```julia
using LyoPronto
```
`LyoPronto` reexports the following packages from `OrdinaryDiffEqRosenbrock`, `OrdinaryDiffEqNonlinearSolve`, `DiffEqCallbacks`, and `Unitful`, so those are available after `using LyoPronto`.

The following packages are not exported by LyoPronto, so you will have to install them and import them as well to follow along here.
```
using DrWatson
using CSV
using TypedTables
using LaTeXStrings
using SavitzkyGolay
using Optim
using Plots
using StatsPlots: @df
using Accessors
```

The following packages are used in this example:

## Loading Experimental Data

First, load a `.csv` of experimental data with `CSV` and `TypedTables`. `CSV.read` method, and use `DrWatson`'s `datadir` function to help keep directory structure clean.

This file has 7 rows of metadata at the top, so the data headers are at row 8; it is located in the project directory under `[project]/data/exp_raw/"2024-06-21-16_MFD_AH.csv`.

```julia
data_raw = CSV.read(datadir("exp_raw", "2024-06-21-16_MFD_AH.csv"), Table, header=8)
```


```julia
pd_data = filter(x->x.phase == 4, procdata)
pd_data.t .-= pd_data.t[1]
```

To check that everything looks right, plot the temperatures, taking advantage of a recipe from this package, as well as the `L"[latex]"` macro from `LaTeXStrings`. We can also exploit the `@df` macro from `StatsPlots` to make this really smooth.

```julia
@df pd_data exptfplot(:t, :T1, :T2, :T3, labels=[L"T_{p1}" L"T_{p2}" L"T_{p3}"])
@df pd_data plot!(:t, :Tsh_sp, c=:black, label=L"T_{sh}")
```

Look at the Pirani pressure and ascertain the end of drying by using a Savitzky-Golay filter to identify a maximum in the second derivative.
Note that because `SavitzkyGolay` doesn't play nice with Unitful, we strip out units and add them back in.
Another plot recipe plots the end of drying as a vertical line on that pressure graph.

```julia
p_pir_sm = savitzky_golay(ustrip.(u"mTorr", pd_data.pch_pir), 91, 3, deriv=0).y *u"mTorr" # 91 a window width; 3 the polynomial order
p_pir_der2 = savitzky_golay(ustrip.(u"mTorr", pd_data.pch_pir), 91, 3, deriv=2).y
t_end = pd_data.t[argmax(p_pir_der2[100:end])+99] # Skip past the beginning, where pressure might be weird
@df pd_data plot(:t, :pch_pir, label="data")
@df pd_data plot!(:t, p_pir_sm, label="smoothed data")
tendplot!(t_end, label=L"t_{end}")
```

Based on an examination of the temperature data, we want to go only up to the "temperature rise" commonly observed in lyophilization near (but not at) the end of drying. This happens at 7.5 hours for `T1` and `T3` and at about 12.5 hours for `T2`.
To pass this information on to the least-squares fitting routine, one could manually trim down all the data to compare, which I have done many times. But since this is such a common operation, I made [a tool](@ref LyoPronto.PrimaryDryFit) which encodes this information:

```julia
fitdat_all = PrimaryDryFit(pd_data.t, 
                        (pd_data.T1[pd_data.t .< 7.5u"hr"],
                        pd_data.T2[pd_data.t .< 12u"hr"],
                        pd_data.T3[pd_data.t .< 7.5u"hr"]),
                        t_end)
```

Note that by passing all three temperature series to `PrimaryDryFit`, this will compare model output to all three temperature series. But in this case, $T_1$ might be an edge vial we want to fit separately, so let's actually use

```julia
fitdat_center = PrimaryDryFit(pd_data.t, 
                        (pd_data.T2[pd_data.t .< 12u"hr"],
                        pd_data.T3[pd_data.t .< 7.5u"hr"]),
                        t_end)
fitdat_edge = PrimaryDryFit(pd_data.t, pd_data.T1[pd_data.t .< 7.5u"hr"])
```
There's also a plot recipe for this type--call `plot(fitdat_center)` and you can get an immediate feel if you got the right data or not.

## Set up model parameters

Below, we make liberal use of Unitful units. Also note that [`RampedVariable`](@ref LyoPronto.RampedVariable) and [`RpFormFit`](@ref LyoPronto.RpFormFit) are used to simplify some common things we use.
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

paramobj = ParamObjPikal(params_bunch)

```

## Run a sanity-check simulation

```julia
# Time span: used to set initial time and to give an upper bound on time, in case parameters are bad
tspan = (0.0, 100.0) # hours
# Initial condition
u0 = ustrip.([u"cm", u"K"], [hf0, Tsh(0u"minute")])

# Set up as an ODE problem
prob = ODEProblem(lyo_1d_dae_f, u0, tspan, paramobj)
# Solve with the Rodas4P() algorithm, use a callback provided by LyoPronto to terminate at end of drying
sol = solve(prob, Rodas4P(), callback=end_drying_callback)

```
And then plot these results with a recipe, to make sure that the chosen $K_{sh-f}$ and $R_p$ are physically reasonable guesses (though they will not be exact):

```julia
@df pd_data plot(:t, :Tsh_sp, c=:black, label=L"T_{sh}")
plot!(fitdat_center)
tendplot!(fitdat_center.t_end, label=L"t_{end}")
modconvtplot!(sol, label=L"$T_p$, model")
```

## Minimize least square difference to compare to experimental data

It's a good idea at this point to make sure that the objective function is working as expected:
```julia
KRp_guess = [12, 0.1, 5, 0.1]
obj_KRp(KRp_guess)
```
should return a number.

Now, to actually carry out the fitting, we find a minimum for this objective function, using `Optim` with the `NelderMead()` algorithm:
```julia
opt_KRp = optimize(obj_KRp, KRp_guess, NelderMead())
```

Get out the found values of our tuning parameters, and generate the corresponding solution profile:

```julia
KRp = Optim.minimizer(opt_KRp)
prof, optparams = gen_sol_conv(KRp)
```

Plot and compare experiment to model:
```julia
@df pd_data exptfplot(:t, :T1, :T2, :T3, labels=[L"T_{p1}" L"T_{p2}" L"T_{p3}"])
plot!(Tsh, label=L"T_{sh}")
modconvtplot!(prof)
tendplot!(t_end)
plot!(xlim = (0, 20))
```

## Going further
Suppose we want to fit a different $K_v$ for the vial I said might be an edge vial, we can write another objective function for that, do new fitting, and add that to the plot.

```julia
# Write a general function that takes K, solves model
function gen_sol_K_edge(K, po::ParamObjPikal, u0, timespan)
    new_params = @set po.Kshf = x->(K*u"W/m^2/K")
    prob = ODEProblem(lyo_1d_dae_f, u0, timespan, new_params)
    sol = solve(prob, Rosenbrock23(autodiff=false), callback=end_drying_callback)
    return sol, new_params
end

# Make a specific version for this experiment, with Rp tuned from above
gen_sol_K = K->gen_sol_K_edge(K, optparams, u0, tspan)

obj_K_edge = Kvec->obj_expT(gen_sol_K(Kvec[1])[1], fitdat_edge)

opt_K_edge = optimize(obj_K_edge, [10.0], NelderMead())

K_edge = Optim.minimizer(opt_K_edge)[1]

prof_edge = gen_sol_K(K_edge)[1]

@df pd_data exptfplot(:t, :T1, :T2, :T3, labels=[L"T_{p1}" L"T_{p2}" L"T_{p3}"])
plot!(Tsh, label=L"T_{sh}", color=:black)
modconvtplot!(prof, prof_edge, labels=[L"$T_{p}$, center" L"$T_{p}$, edge"])
tendplot!(t_end, label=L"t_{end}")
plot!(xlim = (0, 20))
```