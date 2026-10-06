# # Extending Fits with New Experimental Data Types
# The `ExpFitData` type is designed to be extensible, so that other types of experimental data 
# beyond those already implemented can be used in the parameter fitting framework. 
# This is done by defining a new subtype of [`LyoPronto.AbstractExpDatum`](@ref) and 
# implementing the required methods for that type. 
# We will demonstrate this here, with two examples: 
# one which represents a total measurement of ``Q_{vw-f}`` (heat transfer from vial wall to frozen product)
# across all of drying, and one which represents a measurement of this heat transfer at specific time points.
# !!! note Synthetic example
#     These types of data are not necessarily easy to measure or useful to 
#     lyophilization design, but these will work as simple examples to demonstrate
#     the code's interface. 

using LyoPronto
## We will be extending the following functions with new methods, so we need to `import` them
import LyoPronto: AbstractExpDatum, resid_name, obj_exp_datum, 
    err_exp_datum!
## Other packages we need
using TypedTables
using Plots
using StatsPlots: @df


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
@df qs areastackplot(t, :Q_shf, :Q_vwf, :Q_RF_f, c=[:red :orange :yellow], labels=["Shelf" "Vial wall" "Volumetric"], 
    xlabel="Time", ylabel="Heat transfer", title="Heat transfer modes over time", unitformat=:square)

# We will use this base solution to create some synthetic experimental data

# First, do a Riemann integration to get total Q_vwf
t_weights = (vcat(0u"hr", diff(t)) + vcat(diff(t), 0u"hr"))/2
total_Qvwf = sum(qs.Q_vwf .* t_weights)

# Now, evaluate Q_vwf at specific time points
t_exp = range(0.0u"hr", stop = 5.0u"hr", step=10u"minute")
t_exp_nd = ustrip.(u"hr", t_exp)
Q_vwf_exp = [calc_md_Q(base_sol(ti), po, ti).Q_vwf for ti in t_exp_nd]

# ## Set up the experimental data types
# We will define two new types of experimental data, one for the total Q_vwf and one for the time series of Q_vwf.

# The `TotalQvwf` type will represent a single measurement of the total heat transfer from the vial wall to the frozen product over the entire drying process. 
# All it needs to store is a single value of Q_vwf.
struct TotalQvwf{T} <: AbstractExpDatum
    Q_vwf::T
end

# Now, we need to define the required methods for this type.
LyoPronto.time_bound_data(::TotalQvwf) = false
LyoPronto.resid_name(::TotalQvwf) = :Qvwf_total
# Define a residual (for nonlinear least squares fitting) and a loss (for optimization-solver least squares fitting).
## Loss function: just returns the error
function LyoPronto.obj_exp_datum(sol, st, dat::TotalQvwf)
    ## `st` gets computed internally and used to evaluate the solution only at time points
    ## where model and experiment are both available
    q_series = model_result(sol, st; var=:Q_vwf)
    ## Trapezoidal integration over time, really concisely
    tweight = (x->(vcat(0, x)+vcat(x,0))/2)(diff(sol.t))*u"hr"
    qinteg = sum(q_series .* tweight) |> u"W*hr"
    ## The residual is the difference between the model's total Q_vwf and the experimental value
    return (qinteg - dat.Q_vwf)^2
end
## Residual function: fills the error in place in an array, returns number of filled indices
function LyoPronto.err_exp_datum!(err_array, i0, sol, st, dat::TotalQvwf, weight)
    q_series = model_result(sol, st; var=:Q_vwf)
    ## Trapezoidal integration over time, really concisely
    tweight = (x->(vcat(0, x)+vcat(x,0))/2)(diff(sol.t))*u"hr"
    qinteg = sum(q_series .* tweight) |> u"W*hr"
    ## The residual is the difference between the model's total Q_vwf and the experimental value
    ## `io` is the last index that was filled into the array
    err_array[i0+1] = (qinteg - dat.Q_vwf) * weight
    ## The number of residuals filled is 1, since this is a single measurement 
    return 1
end

# The `TimeSeriesQvwf` type will represent a series of measurements of Q_vwf at specific time points.
# Since this may or may not be evenly spaced or at the same time points as other measurements,
# we have fields for a time vector and for a range of indices, but these may be simply set to 
# `missing` if not needed.
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

# Now, we need to define the required methods for this type.
LyoPronto.time_bound_data(::TimeSeriesQvwf) = true
LyoPronto.resid_name(::TimeSeriesQvwf) = :Qvwf
LyoPronto.has_timevec(tq::TimeSeriesQvwf) = !ismissing(tq.t)
LyoPronto.nontrivial_t_range(tq::TimeSeriesQvwf) = tq.t_range != eachindex(tq.Q_vwf)
# Define a residual (for nonlinear least squares fitting) and a loss (for optimization-solver least squares fitting).
## Loss function: just returns the error
function LyoPronto.obj_exp_datum(sol, st, dat::TimeSeriesQvwf)
    ## `st` gets computed internally and used to evaluate the solution only at time points
    ## where model and experiment are both available
    q_series = model_result(sol, st; var=:Q_vwf)
    ## Sum of squared differences between model and experimental
    return sum(abs2, q_series - dat.Q_vwf)
end
## Residual function: fills the error in place in an array, returns number of filled indices
function LyoPronto.err_exp_datum!(err_array, i0, sol, st, dat::TotalQvwf, weight)
    ## The `st` object handles the time alignment between the model solution and the experimental data
    q_series = model_result(sol, st; var=:Q_vwf)
    q_errs = (dat.Q_vwf[st.ti_m_start:st.ti_m_end] .- q_series[begin:st.len])/sqrt(st.len)
    ## Total number of experimental data points is the length of the time series
    nq = length(dat.Q_vwf)
    err_array[i0+st.ti_m_start:i0+st.ti_m_end] .= ustrip.(NoUnits, q_errs * weight)
    ## Fill in the rest of the error array with a sentinel value to indicate that these points are not available
    err_array[i0+1:i0+st.ti_m_start] .= 0.0
    err_array[i0+st.ti_m_end+1:i0+ntvw] .= 0.0
    ## Return the number of residuals filled, which is the length of the time series
    return nq
end

# ## Set up an actual experimental fitting problem

# We will put together two experimental data objects, one for each "measurement" we created above.
# To demonstrate the results of each separately, we make two [`ExpFitData`](@ref) objects,
# but in general any number of experimental data objects can be combined together.
efd1 = ExpFitData(t_exp, TotalQvwf(total_Qvwf))
efd2 = ExpFitData(t_exp, TimeSeriesQvwf(Q_vwf_exp))

####TODO: fill out all the weighting, etc.