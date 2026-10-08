# # Extending Fits with New Experimental Data Types
# The `ExpFitData` type is designed to be extensible, so that other types of experimental data 
# beyond those already implemented can be used in the parameter fitting framework. 
# This is done by defining a new subtype of [`LyoPronto.AbstractExpDatum`](@ref) and 
# implementing the required methods for that type. 
# We will demonstrate this here, with two examples: 
# one which represents a total measurement of ``Q_\mathrm{vw-f}`` (heat transfer from vial wall to frozen product)
# across all of drying, and one which represents a measurement of this heat transfer at specific time points.

# !!! note "Synthetic example"
#     These types of data are not necessarily easy to measure or useful to 
#     lyophilization design, but these will work as simple examples to demonstrate
#     the code's interface. 

using LyoPronto
## We will be extending the following functions with new methods, so we need to `import` them
import LyoPronto: AbstractExpDatum, resid_name, obj_exp_datum, 
    err_exp_datum!
## Other packages we need
using LaTeXStrings
using TypedTables
using Plots
using StatsPlots: @df
using Accessors: setproperties
## For fitting
using TransformVariables
using OptimizationOptimJL
using NonlinearSolve


# ## Demonstration case
# For demonstration, we will use the model for microwave-assisted lyophilization, with conditions 
# like those in [Fitting to Microwave Drying](@ref).


vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
A_p = π*rad_i^2  # cross-sectional area inside the vial
A_v = π*rad_o^2 # vial bottom area
m_v = get_vial_mass(vialsize)
## Formulation and fill
c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16.0u"cm*hr*Torr/g"
A2 = 0.01u"1/cm"
Rp = RpFormFit(R0, A1, A2)
Vfill = 5u"mL"
## Heat transfer
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
K_shf = RpFormFit(KC, KP, KD)
## Geometry
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
## RF fit parameters 
Bf = 5.0e8u"Ω/m^2"
Bvw = 1.1e7u"Ω/m^2"
Kvwf = 22.3u"W/K/m^2"
## Controllable inputs
f_RF = 8u"GHz"
Tsh = RampedVariable(uconvert.(u"K", [-40.0u"°C", 10.0u"°C"]), 0.5u"K/minute")
pch = RampedVariable(100u"mTorr")
## Total of 10W nominal power, multiplied by 0.54 to account for system losses
P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

po = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
))

base_sol = solve(ODEProblem(po), LyoPronto.odealg_chunk3)

# Let's take a look at the heat transfer modes over time.
qs = Table(map(base_sol.u, base_sol.t) do u, t
    calc_md_Q(u, po, t)
end)

## LyoPronto provides a slick plot recipe for this, which stacks up the values over time so you can visually compare areas
t = base_sol.t*u"hr"
@df qs areastackplot(t, :Q_shf, :Q_vwf, :Q_RF_f, c=[:red :orange :yellow], labels=[L"Q_\mathrm{sh-f}" L"Q_\mathrm{vw-f}" L"Q_\mathrm{RF-f}"], 
    xlabel="Time", ylabel="Heat transfer", title="Heat transfer modes over time", unitformat=:square)

# We will use this base solution to create some synthetic experimental data.

# First, do a trapezoidal integration to get total ``Q_\mathrm{vw-f}``:
t_weights = (vcat(0u"hr", diff(t)) + vcat(diff(t), 0u"hr"))/2
total_Qvwf = sum(qs.Q_vwf .* t_weights)

# Now, evaluate ``Q_\mathrm{vw-f}`` at specific time points:
t_exp = range(0.0u"hr", stop = 5.0u"hr", step=10u"minute")
t_exp_nd = ustrip.(u"hr", t_exp)
Q_vwf_exp = [calc_md_Q(base_sol(ti), po, ti).Q_vwf for ti in t_exp_nd]

# ## Set up the experimental data types
# We will define two new types of experimental data, one for the total ``Q_\mathrm{vw-f}`` and one for the time series of ``Q_\mathrm{vw-f}``.

# The `TotalQvwf` type will represent a single measurement of the total heat transfer from the vial wall to the frozen product over the entire drying process. 
# All it needs to store is a single value of ``Q_\mathrm{vw-f}``.
struct TotalQvwf{T} <: AbstractExpDatum
    Q_vwf::T
end

# Now, we need to define the required methods for this type.
LyoPronto.time_bound_data(::TotalQvwf) = false
LyoPronto.resid_name(::TotalQvwf) = :Qvwf_total
LyoPronto.num_errs(::TotalQvwf) = 1
# Define a residual (for nonlinear least squares fitting) and a loss (for optimization-solver least squares fitting).
## Loss function: just returns the error
function LyoPronto.obj_exp_datum(sol, dat::TotalQvwf; verbose=false)
    
    q_series = map(sol.u, sol.t) do u, t
        calc_md_Q(u, sol.prob.p, t).Q_vwf
    end
    ## Trapezoidal integration over time, really concisely
    dt = (x->(vcat(0, x)+vcat(x,0))/2)(diff(sol.t))*u"hr"
    qinteg = sum(q_series .* dt) |> u"W*hr"
    ## Provide some information if requested
    verbose && @info "Total Q_vwf: model = $qinteg, exp = $(dat.Q_vwf)"
    ## The residual is the difference between the model's total Q_vwf and the experimental value
    return (qinteg - dat.Q_vwf)^2
end
## Residual function: fills the error in place in an array, returns number of filled indices
function LyoPronto.err_exp_datum!(err_array, i0, sol, dat::TotalQvwf, weight; verbose=false)
    q_series = map(sol.u, sol.t) do u, t
        calc_md_Q(u, sol.prob.p, t).Q_vwf
    end
    ## Trapezoidal integration over time, really concisely
    tweight = (x->(vcat(0, x)+vcat(x,0))/2)(diff(sol.t))*u"hr"
    qinteg = sum(q_series .* tweight) |> u"W*hr"
    ## Provide some information if requested
    verbose && @info "Total Q_vwf: model = $qinteg, exp = $(dat.Q_vwf)"
    ## The residual is the difference between the model's total Q_vwf and the experimental value
    ## `io` is the last index that was filled into the array
    err_array[i0+1] = ustrip(NoUnits, (qinteg - dat.Q_vwf) * weight)
    ## The number of residuals filled is 1, since this is a single measurement 
    return 1
end

# !!! note "`verbose` keyword argument"
#     The `verbose` keyword argument to `obj_exp_datum` and `err_exp_datum!` 
#     is provided to the loss and residual functions so that,
#     if you are trying to debug a fit,
#     you can see the model's predictions and the experimental data at each iteration of the optimization.
#     This is especially useful if you are trying to figure out why a fit is not converging 
#     or is converging to a poor solution.
#     
#     Make sure you provide it as an allowed keyword argument, because otherwise you will 
#     get an error.
#     Set the default of `verbose=false` so that you don't get a lot of output during normal fitting,
#     but you can set it to `true` when you want to see the output.

# The `TimeSeriesQvwf` type will represent a series of measurements of ``Q_\mathrm{vw-f}`` at specific time points.
# Since this may or may not be evenly spaced or at the same time points as other measurements,
# we have fields for a time vector and for a range of indices, but these may be simply set to 
# `missing` or `1:length(Q_vwf)` respectively if not needed.
# For this example, we will not set those values.
struct TimeSeriesQvwf{T1, T2, T3} <: AbstractExpDatum
    Q_vwf::T1
    t_range::T2
    t::T3
end
## Set up a constructor which fills in default values when t is in common with other measurements
function TimeSeriesQvwf(Q_vwf, t_range=eachindex(Q_vwf); t=missing)
    ## We could do extra data validation right here
    return TimeSeriesQvwf(Q_vwf, t_range, t)
end

# Now, we need to define the required interface methods for this type.
LyoPronto.time_bound_data(::TimeSeriesQvwf) = true # this data is bound to specific time points
LyoPronto.resid_name(::TimeSeriesQvwf) = :Qvwf # provide a symbol which will map to weighting in residual
LyoPronto.num_errs(q::TimeSeriesQvwf) = length(q.Q_vwf) # indicate how many residuals the data can have
LyoPronto.has_timevec(tq::TimeSeriesQvwf) = !ismissing(tq.t) # indicates whether a separate time vector is provided
LyoPronto.nontrivial_t_range(tq::TimeSeriesQvwf) = tq.t_range != eachindex(tq.Q_vwf) # indicates whether the time range is a nontrivial subset of the time vector



# ### Define a loss function (for optimization-solver least squares fitting)
# Since this is a time series, we will need to interpolate the model solution to the time points of the experimental data.
# This is handled internally and managed by the `st` object, which the `obj_exp` function will create and pass to this function.
# You will need to call `LyoPronto.exp_time_inds(st)` to get the experimental time indices,
# `LyoPronto.model_time_inds(st)` to get the model time indices,
# and `length(st)` to get the length of the time series.
# These will be used to index into the experimental data and the model solution, respectively, 
# and the length is used to normalize the residuals so that the residuals are approximately 
# independent of the number of time points which could actually be used in the fit.
## Loss function: just returns the error
function LyoPronto.obj_exp_datum(sol, st, dat::TimeSeriesQvwf; verbose=false)
    ## `st` gets computed internally and used to evaluate the solution only at time points
    ## where model and experiment are both available
    q_series = model_result(sol, st; var=:Q_vwf)
    ## Provide some information if requested
    verbose && @info "Q_vwf:" sum(q_series) LyoPronto.exp_time_inds(st)
    ## Mean sum of squared differences between model and experimental
    return sum(abs2, q_series[LyoPronto.model_time_inds(st)] - dat.Q_vwf[LyoPronto.exp_time_inds(st)])/length(st)
end
# ### Define a residual function (for nonlinear least squares fitting)
## Residual function: fills the error in place in an array, returns number of filled indices
function LyoPronto.err_exp_datum!(err_array, i0, sol, st, dat::TimeSeriesQvwf, weight; verbose=false)
    ## The `st` object handles the time alignment between the model solution and the experimental data
    q_series = model_result(sol, st; var=:Q_vwf)
    ## Compute the errors, and normalize by the `sqrt(size of time series)`
    ## (`sqrt` because if these are squared and summed, it will correspond to above loss)
    q_errs = (dat.Q_vwf[LyoPronto.exp_time_inds(st)] .- q_series[LyoPronto.model_time_inds(st)])/sqrt(length(st))
    ## Total number of experimental data points is the length of the time series
    nq = length(dat.Q_vwf)
    ## Provide some debug info if requested
    verbose && @info "Q_vwf errors:" q_errs LyoPronto.exp_time_inds(st)
    ## Fill in the portion of the error array which corresponds to these data 
    ## with a sentinel value (which indicates the model couldn't evaluate)
    err_array[i0 .+ (1:nq)] .= 0.0
    ## then overwrite the portion where model values are available with the actual errors, weighted appropriately
    err_array[i0 .+ LyoPronto.exp_time_inds(st)] .= ustrip.(NoUnits, q_errs * weight)
    ## Return the number of residuals filled, which is the length of the time series
    return nq
end


# ## Weighting

# !!! tip "Dimensional weighting"
#     Weighting in the loss function is unavoidable if we provide experimental data with different dimensions 
#     (e.g., because we can't add watts to kelvins). 
#     LyoPronto uses the `Unitful` package to enforce dimensional consistency, 
#     so if you neglect to add dimensional weighting for your loss function, 
#     you will get an error.
# 
#     In practice you will usually define it together with your fitting problem,
#     since you may need to tinker 
#     with the weights to get the most useful fit:
#     for example, do you care more about the model fitting the temperatures or the total drying time?

# For a loss function, 
# the dimensions of the weight need to be the inverse square of the residual dimension,
# so that the weighted loss is dimensionless.
# For a residual function, 
# the dimensions of the weight need to be the inverse of the residual dimension,
# so that the weighted residual is dimensionless.

# The keyword argument name to [`loss_weighting`](@ref)
# or to [`residual_weighting`](@ref) should match the symbol 
# returned by `resid_name` for the data type (as defined above); 
# the default keyword arguments (`t`, `Tf`, and `Tvw`) correspond to `resid_name(::EndTimeData)`,
# `resid_name(::TfData)`, and `resid_name(::TvwSeriesData)` or `resid_name(::TvwEndData)`, 
# respectively, so you should also decide appropriate weights for those keyword arguments 
# if you are using those data types in your fitting problem.

# To demonstrate, we will construct some example weightings for the two experimental data types we just defined.
single_q_loss_weights = loss_weighting(Qvwf_total = 1.0u"W^-2*hr^-2", Tvw=10.0u"K^-2")
multi_q_loss_weights = loss_weighting(Qvwf = 1.0/u"W^2")
single_q_resid_weights = residual_weighting(Qvwf_total = 1.0/u"W*hr", Tf=0.5u"K^-1")
multi_q_resid_weights = residual_weighting(Qvwf = 1.0*u"W^-1", t=0.01u"hr^-1")

# ## Set up the fitting problem

# We will put together two experimental data objects, one for each "measurement" we created above.
# To demonstrate the results of each separately, we make two [`ExpFitData`](@ref) objects,
# but in general any number of experimental data objects can be combined together.
efd1 = ExpFitData(t_exp, TotalQvwf(total_Qvwf))
efd2 = ExpFitData(t_exp, TimeSeriesQvwf(Q_vwf_exp))

# As a toy problem, we will say that we are confident in our estimate of the parameter 
# ``K_\mathrm{vw-f}``, 
# but not ``B_\mathrm{f}`` and ``B_\mathrm{vw}``,
# so we will adjust ``B_\mathrm{f}`` and ``B_\mathrm{vw}`` to fit the synthetic experimental data.

# We map from dimensionless parameter vector to the named physical parameters
# using `TransformVariables` so that the optimization solver can work in an unbounded dimensionless space.
trans = as((; 
    Bf = TVScale(1e9u"Ω/m^2") ∘ TVLogistic(), #Bracket from 0 to 1e9 Ω/m^2
    Bvw = TVScale(1e8u"Ω/m^2") ∘ TVLogistic(), #Bracket from 0 to 1e8 Ω/m^2
))


# ## Optimization problem

# We will use the `OptimizationOptimJL` package to solve the optimization problem,
# which is a wrapper around the `Optim` package.

# Our optimization functions are as follows, using [`obj_pd`](@ref) to compute the objective function
# and the `ForwardDiff` package for automatic differentiation of the objective function.
obj_single_q = OptimizationFunction((u, p) -> obj_pd(u, p; weights=single_q_loss_weights), AutoForwardDiff())
obj_multi_q = OptimizationFunction((u, p) -> obj_pd(u, p; weights=multi_q_loss_weights), AutoForwardDiff())

p0 = [-3.0, -3.0] # initial guess in dimensionless space
transform(trans, p0) # check the values in dimensional space of our initial guess

# It's a good idea at this point to make sure that initial guesses are reasonable.
# For the microwave model, one quantitative way is to check that the initial guesses
# satisfy an energy balance on electromagnetic terms:
guessed_po = setproperties(po, transform(trans, p0))
LyoPronto.rf_lumcap_EM_violate(guessed_po) # should be false to be reasonable

# If we had used `[0.0, 0.0]` as an initial guess, the energy balance would be violated: 
guessed_po = setproperties(po, transform(trans, [0.0, 0.0]))
LyoPronto.rf_lumcap_EM_violate(guessed_po) # should be false to be reasonable

# The initial guess gets passed to the optimization problem, 
# along with the transformation, other experimental conditions in `po`, and the experimental data.
# Then, we solve the optimization problem using the `BFGS` algorithm from `Optim`.
opt1 = solve(OptimizationProblem(obj_single_q, p0, (trans, po, efd1)), Optim.BFGS())
# For the time series data, we use the same initial guess and transformation, but with the second experimental data object.
opt2 = solve(OptimizationProblem(obj_multi_q, p0, (trans, po, efd2)), Optim.BFGS())

# Let's see how close they got to matching the total ``Q_\mathrm{vw-f}``:
(opt1.objective, opt2.objective) # should be close to zero if the fit was good
# Both very close. 

# ## Nonlinear least squares problem
# We can also set up a nonlinear least squares problem using the `NonlinearSolve` package.
# There are some theoretical advantages to this, since the optimization problem has to take
# an extra derivative to compute essentially these residuals, 
# but the difficulty for us is in handling the variable number of residuals available 
# from the model
# (which depends on the model parameters, because those will change total drying time).

# Set up a nonlinear least squares problem using a similar approach to the optimization problem, 
# but [`nls_pd!`](@ref) will be used to compute the residuals; 
# since the length of a cached residual array needs to match the number of experimental data points,
# LyoPronto provides a convenience constructor for the `NonlinearFunction` with the signature 
# `NonlinearSolve.NonlinearFunction(::ExpFitData; kwargs...)`.

# Here, we will pass a keyword argument `badprms` to enforce 
# the electromagnetic energy balance we used above to check if our initial guess was reasonable.
# This can be done for the optimization problem as well, by passing `badprms=...` in the same place as `weights=...`.
nlf1 = NonlinearFunction(efd1; weights=single_q_resid_weights, badprms=LyoPronto.rf_lumcap_EM_violate)
nlf2 = NonlinearFunction(efd2; weights=multi_q_resid_weights, badprms=LyoPronto.rf_lumcap_EM_violate)

nls1 = solve(NonlinearLeastSquaresProblem(nlf1, p0, (trans, po, efd1)), LevenbergMarquardt())
# And again for the time series data:
nls2 = solve(NonlinearLeastSquaresProblem(nlf2, p0, (trans, po, efd2)), LevenbergMarquardt())

# ## Compare fit results to original data

# For all these cases, we can compare the original parameters to the fitted values:
casenames = ["Opt, integrated", "Opt, time", "NLS, integrated", "NLS, time"]
cases = [opt1, opt2, nls1, nls2]
markers = [:utriangle :dtriangle :ltriangle :rtriangle]
guess = transform(trans, p0)

plot(xlabel=L"B_\mathrm{f}", ylabel=L"B_\mathrm{vw}", title="Fit results in parameter space", legend=:topleft)
scatter!([guess.Bf], [guess.Bvw], label="Initial guess", ms=10, shape=:square)
for (name, opt, marker) in zip(casenames, cases, markers)
    fit = transform(trans, opt.u)
    scatter!([fit.Bf], [fit.Bvw], label=name, ms=10, shape=marker)
end
scatter!([Bf], [Bvw], label="\"True\"", ms=8, shape=:star5)

# So the fits which had time series data were able to match the true parameters closely
# (indicated by overlapping markers at the "true" values in the plot above),
# while the fits which only had total ``Q_\mathrm{vw-f}`` did not. The NLS solver did not 
# even move from the initial guess, which suggests that some tolerances might need tweaking, but 
# of course the more productive approach is to give better data.

# Next, let's compare the values over time of ``Q_\mathrm{vw-f}`` predicted by the model to the synthetic experimental data we created.

t_qvwfs = map(cases) do opt
    sol = gen_sol_pd(opt.u, trans, po)
    qvwf = map(sol.u, sol.t) do u, t
        calc_md_Q(u, po, t).Q_vwf
    end
    return sol.t*u"hr", qvwf
end

plot(u"hr", u"W", xlabel="Time", ylabel=L"Q_\mathrm{vw-f}", title="Model predictions vs synthetic experimental data")
for (label, (t, qvwf), marker) in zip(casenames, t_qvwfs, markers)
    plot!(t, qvwf; label, marker, lw=2, fillto=0, fillalpha=0.1)
    @info "check" name t qvwf
end
scatter!(t_exp, Q_vwf_exp, label="Synthetic experimental data", color=:black)

# So, both the optimization and nonlinear least squares did a good job of fitting
# when provided the full time series of ``Q_\mathrm{vw-f}``, 
# but unsurprisingly, when the only provided data was matching the total ``Q_\mathrm{vw-f}``,
# the fit doesn't look as good in time, nor does it match the true parameters.

# ## Adding end time to total Qvwf to improve fit

# So, as a scientific curiosity, can we get a good overall fit if we match 
# total ``Q_\mathrm{vw-f}`` and also match the drying time?

t_end = base_sol.t[end]*u"hr"
efd3 = ExpFitData(t_exp, TotalQvwf(total_Qvwf), EndTimeData(t_end))

# Since `single_q_loss_weights` already has a weight for `t`, we don't need a new objective function.
# It's enough to add the new experimental data.

opt3 = solve(OptimizationProblem(obj_single_q, p0, (trans, po, efd3)), Optim.BFGS())

(obj=opt3.objective, Bf=transform(trans, opt3.u).Bf, Bvw=transform(trans, opt3.u).Bvw)

# Let's see how close it got to matching the time series of ``Q_\mathrm{vw-f}`` and the drying time:

sol3 = gen_sol_pd(opt3.u, trans, po)
qvwf3 = map(sol3.u, sol3.t) do u, t
    calc_md_Q(u, po, t).Q_vwf
end
t3 = sol3.t*u"hr"

plot(u"hr", u"W", xlabel="Time", ylabel=L"Q_\mathrm{vw-f}", title="Model predictions vs synthetic experimental data")
plot!(t3, qvwf3; label=L"Both $t_\mathrm{end}$ and total $Q_\mathrm{vw-f}$", marker=:diamond, lw=2, fillto=0, fillalpha=0.1)
scatter!(t_exp, Q_vwf_exp, label="Synthetic experimental data", color=:black)

# So providing the end time as an additional experimental data point 
# improved things significantly, even though it didn't yield an exact match.
# This principle will hold in general: the more physical data we can incorporate into the 
# fitting problem, the better our results will be.
# To put this in machine learning terms, the more features we can provide to the model, 
# the better it will be able to approximate the true physics.
