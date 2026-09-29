
# Data objects for ExpFitData. Each holds one kind of experimental measurement
# plus an optional per-object time vector (`t`). When `t` is `missing`, the
# container's `t` is used (see `fit_t`).
abstract type AbstractExpDatum end

# Trait for identifying whether a given fit type has an associated vector of times
function has_timevec end
# Trait for identifying dimensions of the residuals for a given fit type
function error_dims end

"""
    TfData(Tf; t=missing)
    TfData(Tf, t_range; t=missing)

A single experimental freezing-front (Tf) temperature series, for use in [`ExpFitData`](@ref).

Fields:
- `Tf`: an `AbstractVector` of `Temperature` (one series).
- `t_range`: a range of integer indices into the time vector that `Tf` corresponds to;
  `Tf[k]` is the temperature measured at time `t[t_range[k]]`.
- `t`: an optional `AbstractVector` of `Time`; when `missing` (the default), the container's
  time vector is used (see [`fit_t`](@ref)).

Constructors:
- `TfData(Tf; t=missing)`: the time-index range defaults to `1:length(Tf)`, i.e. the series
  corresponds to the first `length(Tf)` points of the time vector.
- `TfData(Tf, t_range; t=missing)`: use an explicit `t_range` (a range of the same length
  as `Tf`), e.g. to indicate the series was measured over a sub-window of the time vector.
"""
struct TfData{T1, T2, T3} <: AbstractExpDatum
    Tf::T1        # AbstractVector of Temperature (one series)
    t_range::T2  # range of time indices that Tf corresponds to
    t::T3          # AbstractVector of Time; `missing` → use container `t`
    function TfData(Tf, t_range; t=missing)
        if !(Tf isa AbstractVector) || !(first(Tf) isa Unitful.Temperature)
            throw(ArgumentError("Tf should be a vector of temperature measurements"))
        end
        if !(t_range isa AbstractRange) || !(first(t_range) isa Integer)
            throw(ArgumentError("t_range should be a range of integer time indices"))
        end
        if length(t_range) != length(Tf)
            throw(ArgumentError("t_range must have the same length as Tf"))
        end
        if !ismissing(t) && (!(t isa AbstractVector) || !(first(t) isa Unitful.Time))
            throw(ArgumentError("t should be either `missing` or a vector of time points with units of time"))
        end
        new{typeof(Tf), typeof(t_range), typeof(t)}(Tf, t_range, t)
    end
end

# Convenience constructor: default the time-index range to the full series length.
function TfData(Tf; t=missing)
    return TfData(Tf, 1:length(Tf); t)
end

has_timevec(::TfData) = true

# Provide separate types for Tvw as a full series or as a single endpoint
"""
    TvwSeriesData(Tvws; t=missing)
    TvwSeriesData(Tvws, t_range; t=missing)

A single experimental vial-wall (Tvw) temperature series, for use in [`ExpFitData`](@ref).

Fields:
- `Tvws`: an `AbstractVector` of `Temperature` (one series).
- `t_range`: a range of integer indices into the time vector that `Tvws` corresponds to;
  `Tvws[k]` is the temperature measured at time `t[t_range[k]]`.
- `t`: an optional `AbstractVector` of `Time`; when `missing` (the default), the container's
  time vector is used (see [`fit_t`](@ref)).

Constructors:
- `TvwSeriesData(Tvws; t=missing)`: the time-index range defaults to `1:length(Tvws)`, i.e. the
  series corresponds to the first `length(Tvws)` points of the time vector.
- `TvwSeriesData(Tvws, t_range; t=missing)`: use an explicit `t_range` (a range of the
  same length as `Tvws`), e.g. to indicate the series was measured over a sub-window of the
  time vector.
"""
struct TvwSeriesData{T1, T2, T3} <: AbstractExpDatum
    Tvw::T1       # AbstractVector of Temperature
    t_range::T2 # range of time indices that Tvws corresponds to
    t::T3          # AbstractVector of Time; `missing` → use container `t`
    function TvwSeriesData(Tvw, t_range; t=missing)
        if !(Tvw isa AbstractVector) || !(first(Tvw) isa Unitful.Temperature)
            throw(ArgumentError("Tvw should be a temperature or a vector of temperature measurements"))
        end
        if !(t_range isa AbstractRange) || !(first(t_range) isa Integer)
            throw(ArgumentError("t_range should be a range of integer time indices"))
        end
        if length(t_range) != length(Tvw)
            throw(ArgumentError("t_range must have the same length as Tvw"))
        end
        if !ismissing(t) && (!(t isa AbstractVector) || !(first(t) isa Unitful.Time))
            throw(ArgumentError("t should be a vector of time points with units of time"))
        end
        new{typeof(Tvw), typeof(t_range), typeof(t)}(Tvw, t_range, t)
    end
end

# Convenience constructor: default the time-index range to the full series length.
function TvwSeriesData(Tvw; t=missing)
    return TvwSeriesData(Tvw, 1:length(Tvw); t)
end

has_timevec(::TvwSeriesData) = true

"""
    TvwEndData(Tvw_end)
A single experimental vial-wall (Tvw) temperature endpoint, for use in [`ExpFitData`](@ref).
"""
struct TvwEndData{T1} <: AbstractExpDatum
    Tvw_end::T1     # Unitful.Temperature
    function TvwEndData(Tvw_end)
        if !(Tvw_end isa Unitful.Temperature)
            throw(ArgumentError("Tvw_end should be a temperature"))
        end
        new{typeof(Tvw_end)}(Tvw_end)
    end
end

has_timevec(::TvwEndData) = false

"""
    EndTimeData(t_end)
An experimentally measured drying time, for use in [`ExpFitData`](@ref). 

`t_end` can be a single time or a tuple of two times, indicating a window of acceptable drying times.
"""
struct EndTimeData{T1} <: AbstractExpDatum
    t_end::T1      # Unitful.Time OR Tuple of two Times
    function EndTimeData(t_end)
        if !(t_end isa Unitful.Time) && !(t_end isa Tuple && first(t_end) isa Unitful.Time)
            throw(ArgumentError("t_end should be a time or a tuple of two times"))
        end
        new{typeof(t_end)}(t_end)
    end
end

has_timevec(::EndTimeData) = false

# --------
# Collector struct for all data to be fit, in a single experiment

struct ExpFitData{T1, T2}
    t::T1
    data::T2       # Tuple of TfData, TvwSeriesData, TvwEndData, EndTimeData objects
    function ExpFitData(t, data::Tuple)
        if !(t isa AbstractVector) || !(first(t) isa Unitful.Time)
            throw(ArgumentError("t should be a vector of time points with units of time"))
        end
        if isempty(data)
            throw(ArgumentError("data should be a non-empty tuple of data objects"))
        end
        for obj in data
            if !(obj isa AbstractExpDatum)
                throw(ArgumentError("data should be a tuple of TfData, TvwSeriesData, TvwEndData, or EndTimeData objects"))
            end
        end
        new{typeof(t), typeof(data)}(t, data)
    end
end

# Order-agnostic convenience constructor: accepts one or more pre-built data
# objects in any order. The order of the arguments determines the order of the
# `data` tuple (and hence the order of residuals in the fit).
function ExpFitData(t, objs::AbstractExpDatum...)
    if isempty(objs)
        throw(ArgumentError("at least one data object is required"))
    end
    return ExpFitData(t, objs)
end

"""
    fit_t(fitdat, obj)

Return the time vector to use for `obj`: `obj.t` if it is set, otherwise the
container's `fitdat.t`.
"""
function fit_t(fitdat::ExpFitData, obj::AbstractExpDatum)
    if has_timevec(obj) 
        return ismissing(obj.t) ? fitdat.t[obj.t_range] : obj.t[obj.t_range]
    else
        return fitdat.t
    end
end

@doc """
    ExpFitData(t, data)   

A type for indicating how experimental data should be fit.

`ExpFitData` is a container holding a master time vector `t` and a tuple `data`
of experimental data objects: any number of [`TfData`](@ref), [`TvwSeriesData`](@ref),
[`TvwEndData`](@ref), and [`EndTimeData`](@ref) objects. The order of objects in `data`
determines the order of residuals in the fitting procedure.

Provided constructors:

    ExpFitData(t, data)          # data is a tuple of pre-built data objects
    ExpFitData(t, obj1, obj2, ...)  # convenience: one or more pre-built data objects, in any order

The convenience constructors that accept raw temperature vectors with keyword
arguments (`Tvws=`, `t_end=`) are only available through the deprecated
[`PrimaryDryFit`](@ref) function.

Each `TfData` and `TvwSeriesData` object has a `t_range` field: the range of time indices into
the time vector that the series corresponds to, so that `series[k]` is the temperature measured
at `t[t_range[k]]`. By default `t_range` is `1:length(series)`, but it can be set to a sub-window
to indicate the series was measured over only part of the time vector. A single endpoint vial wall
temperature is stored in a `TvwEndData` object (no `t_range` needed).

`t_end` (via a `EndTimeData` object) indicates an end of drying, as would be determined from
non-temperature measurements (e.g. Pirani-CM convergence). If no `EndTimeData` object is provided, it is
ignored in the objective function. If set to a tuple of two times, then in the objective function
any time in that window is not penalized; outside that window, squared error takes over, as for
the single time case.

Each data object may carry its own `t` (time vector). When `t` is `missing` (the default in the
data object constructors), the container's `t` is used. This allows mixing data measured on
different time grids.

Common Cases:
- Conventional, single thermocouple: `ExpFitData(t, TfData(Tf))`
- Conventional, multiple thermocouples: `ExpFitData(t, (TfData(Tf1), TfData(Tf2), ...))`
- Conventional with Pirani ending: `ExpFitData(t, (TfData(Tf), EndTimeData(t_end)))`
- RF with measured vial wall: `ExpFitData(t, (TfData(Tf), TvwSeriesData(Tvws)))`
- RF, matching model Tvw to experimental Tf[end] without measured vial wall: `ExpFitData(t, (TfData(Tf), TvwEndData(Tvw_end)))`
"""
ExpFitData

# Convenience constructor for the old PrimaryDryFit API.
# These build TfData/TvwSeriesData/TvwEndData/EndTimeData objects with `t = missing` (use container t).
"""
    $(SIGNATURES)

Deprecated: this is now a constructor for [`ExpFitData`](@ref).
"""
function PrimaryDryFit(t, Tfs; Tvws=missing, t_end=missing)
    Base.depwarn("PrimaryDryFit is deprecated; use ExpFitData instead. The PrimaryDryFit API is now a thin alias over ExpFitData, and will be removed in future.", :PrimaryDryFit)
    # Tfs is a vector or tuple of vectors; build data objects
    Tfs = Tfs
    if Tfs isa AbstractVector
        if eltype(Tfs) <: Number
            Tfs = (Tfs,)
        else
            Tfs = Tuple(Tfs...)
        end
    end
    if Tvws isa AbstractVector
        if eltype(Tvws) <: Number
            Tvws = (Tvws,)
        else
            Tvws = Tuple(Tvws...)
        end
    end
    if t_end isa Tuple
        t_end = extrema(t_end)
    end
    objs = Tuple(TfData(Tf) for Tf in Tfs)
    if !ismissing(Tvws)
        tvw_obj = Tvws isa Number ? TvwEndData(Tvws) : Tuple(TvwSeriesData(Tvw) for Tvw in Tvws)
        objs = (objs..., tvw_obj...)
    end
    if !ismissing(t_end)
        objs = (objs..., EndTimeData(t_end))
    end
    data = objs
    ExpFitData(t, data)  # calls the actual current constructor
end
PrimaryDryFit(t, Tfs, Tvws, t_end) = PrimaryDryFit(t, Tfs; Tvws=Tvws, t_end=t_end)


# TODO: do this programmatically for all AbstractExpDatum types, rather than hard-coding each one
function Base.:(==)(a::TfData, b::TfData)
    return a.Tf == b.Tf && a.t_range == b.t_range && (ismissing(a.t) == ismissing(b.t))
end
function Base.:(==)(a::TvwSeriesData, b::TvwSeriesData)
    return a.Tvw == b.Tvw && a.t_range == b.t_range && (ismissing(a.t) == ismissing(b.t))
end
function Base.:(==)(a::TvwEndData, b::TvwEndData)
    return a.Tvw_end == b.Tvw_end
end
function Base.:(==)(a::EndTimeData, b::EndTimeData)
    return a.t_end == b.t_end
end
function Base.:(==)(a::ExpFitData, b::ExpFitData)
    return a.t == b.t && a.data == b.data
end

function Base.show(io::IO, d::TfData)
    print(io, "TfData(", length(d.Tf), " pts, t_range=", d.t_range, ")")
end
function Base.show(io::IO, d::TvwEndData)
    print(io, "TvwEndData(", d.Tvw_end, ")")
end
function Base.show(io::IO, d::TvwSeriesData)
    print(io, "TvwSeriesData(", length(d.Tvw), " pts, t_range=", d.t_range, ")")
end
function Base.show(io::IO, d::EndTimeData)
    print(io, "EndTimeData(", d.t_end, ")")
end
function Base.show(io::IO, d::ExpFitData)
    # compact = get(io, :compact, false)::Bool
    # if compact
    #     return print(io, "ExpFitData(...)")
    # end
    kinds = [string(typeof(obj).name.name) for obj in d.data]
    print(io, "ExpFitData(", length(d.t), " t pts; ", join(kinds, ", "), ")")
end

# --- Shared solution-trimming helpers -------------------------------------
# These extract the model temperature series (index `idx`) on the experimental
# time grid `t`, handling the pre-interpolated vs. interpolated cases and the
# sub-zero interpolation correction. `i_solstart` is the index into `t` where
# the model solution begins.

struct SolTrim{T}
    i_solstart::Int
    preinterp::Bool
    tmd::T
end

"""
    trim_sol(sol::ODESolution, t_exp)

Trim the solution `sol` to the experimental time grid `t_exp`, returning a `SolTrim` object.

This should be called separately for data objects which have their own time vector.
"""
function trim_sol(sol::ODESolution, t_exp)
    tmd = sol.t[end].*u"hr"
    nt = length(sol.t) - 1
    i_solstart = searchsortedfirst(t_exp, sol.t[begin]*u"hr")
    # Identify if the solution is pre-interpolated to the time points in t
    preinterp = mapreduce(≈, &, sol.t[1:nt].*u"hr", t_exp[i_solstart:i_solstart+nt-1]) 
    SolTrim(i_solstart, preinterp, tmd)
end

function model_result(sol::ODESolution, st::SolTrim, idx; verbose=false)
    if st.preinterp
        return sol[idx, begin:end-1].*u"K" # Leave off last time point because is end time
    end
    tmd = st.tmd
    trim = sol.t[begin]*u"hr" .< st.t .< tmd
    t_trim = st.t[trim]
    Tmd = sol.(ustrip.(u"hr", t_trim), idxs=idx).*u"K"
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tmd .< 0u"K")
        subzero = findall(Vector(Tmd .< 0u"K"))
        Tmd[subzero] .= Tmd[subzero[1] - 1]
        verbose && @info "bad interpolation" subzero Tmd[subzero]
    end
    return Tmd
end

# --- Per-data-type objective functions ------------------------------------
# Each returns a scalar contribution (in K^2 for temperature, hr^2 for time).

function obj_Tf(sol::ODESolution, obj::TfData, t; verbose=false)
    st = trim_sol(sol, t)
    obj_Tf(sol, st, obj; verbose)
end
function obj_Tf(sol::ODESolution, st::SolTrim, dat::TfData; verbose=false)
    Tmd = model_result(sol, st, 2; verbose) # Tf at index 2
    maxind = min(length(dat.t_range), length(Tmd))
    return sum(abs2, (dat.Tf[st.i_solstart:maxind] .- Tmd[begin:maxind-st.i_solstart+1]))/(maxind-st.i_solstart+1)
end

function obj_Tvw(sol::ODESolution, obj::TvwSeriesData, t; verbose=false)
    st = trim_sol(sol, t)
    obj_Tvw(sol, st, obj; verbose)
end
function obj_Tvw(sol::ODESolution, st::SolTrim, dat::TvwSeriesData; verbose=false)
    Tmd = model_result(sol, st, 3; verbose) # Tvw at index 3
    maxind = min(length(dat.t_range), length(Tmd))
    return sum(abs2, (dat.Tvw[st.i_solstart:maxind] .- Tmd[begin:maxind-st.i_solstart+1]))/(maxind-st.i_solstart+1)
end

function obj_Tvw(sol::ODESolution, obj::TvwEndData; verbose=false)
    return (sol[3, end]*u"K" - uconvert(u"K", obj.Tvw_end))^2
end

function obj_tend(sol::ODESolution, obj::EndTimeData; verbose=false)
    tmd = sol.t[end].*u"hr"
    t_end = obj.t_end
    if t_end isa Tuple # See if is inside window and scale appropriately
        mid_t = (t_end[1] + t_end[2]) / 2.0
        if tmd < t_end[1]
            return (mid_t - tmd)^2
        elseif tmd > t_end[2]
            return (mid_t - tmd)^2
        else # Inside window, so no error
            return 0.0u"hr^2"
        end
    else # Compare to a single drying time
        return (t_end - tmd)^2
    end
end

# Function which dispatches to appropriate per-data-type objectives
obj_exp_datum(sol::ODESolution, st::SolTrim, obj::TfData; verbose=false) = obj_Tf(sol, st, obj; verbose)
obj_exp_datum(sol::ODESolution, st::SolTrim, obj::TvwSeriesData; verbose=false) = obj_Tvw(sol, st, obj; verbose)
obj_exp_datum(sol::ODESolution, obj::TfData, t; verbose=false) = obj_Tf(sol, obj, t; verbose)
obj_exp_datum(sol::ODESolution, obj::TvwSeriesData, t; verbose=false) = obj_Tvw(sol, obj, t; verbose)
obj_exp_datum(sol::ODESolution, obj::TvwEndData; verbose=false) = obj_Tvw(sol, obj; verbose)
obj_exp_datum(sol::ODESolution, obj::EndTimeData; verbose=false) = obj_tend(sol, obj; verbose)
# Catch the case where solution failed
obj_exp_datum(sol::Val{NaN}, obj::AbstractExpDatum; verbose=false) = Inf
obj_exp_datum(sol::Val{NaN}, obj::AbstractExpDatum, t; verbose=false) = Inf
obj_exp_datum(sol::Val{NaN}, st::SolTrim, obj::AbstractExpDatum; verbose=false) = Inf

"""
    $(SIGNATURES)

Evaluate an objective function which compares model solution computed by `sol` to experimental data in `efd`.

- `sol` is a solution to an appropriate model; see [`gen_sol_pd`](@ref) for a helper.
- `efd` is an instance of [`ExpFitData`](@ref), which contains some information about what to compare.
- `tweight = 1.0u"K^2/hr^2"` gives the weighting (should have dimensions like K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
- `Tvw_weight = 1.0` gives the weighting of Tvw in the objective, as compared to Tf.

Note that if `efd` has vial wall temperatures (i.e. a [`TvwSeriesData`](@ref) or [`TvwEndData`](@ref) object in `efd.data`), the third-index variable in `sol` is assumed to be temperature, as is true for the lumped capacitance model (see [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `efd`, squared error is computed for each separately then summed; likewise for `Tvw`.
"""
function obj_exp(sol::ODESolution, efd::ExpFitData;
    tweight=1.0u"K^2/hr^2", verbose = false, Tvw_weight=1.0)
    if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
        verbose && @warn "ODE solve did not reach end of drying. Either parameters are bad, or tspan is not large enough." sol.retcode sol.prob.p.hf0 sol[end]
        return Inf
    end
    Tfobj = 0.0u"K^2"
    Tvwobj = 0.0u"K^2"
    tobj = 0.0u"hr^2"
    # Check which, if any, of the data objects have their own time vector provided
    separate_time = map((x->has_timevec(x) && !ismissing(x.t)), efd.data)
    container_trim = trim_sol(sol, efd.t)
    # Iterate over all data objects, providing appropriate solution trimming if necessary
    # It would be better if this were a map of some sort, but splitting up residual types
    # makes that difficult
    for (obj, sep_trim) in zip(efd.data, separate_time)
        st = if sep_trim
            trim_sol(sol, fit_t(efd, obj))
        elseif has_timevec(obj)
            container_trim
        else
            nothing
        end
        # This is ugly. I would like to do this with dispatch, but that probably means some 
        # sort of weighting object that carries the weights around and new types for each 
        # residual or for all residuals together, which feels annoying, so here we are.
        if obj isa TfData
            Tfobj += obj_exp_datum(sol, st, obj; verbose=verbose)
        elseif obj isa TvwSeriesData 
            Tvwobj += obj_exp_datum(sol, st, obj; verbose=verbose)
        elseif obj isa TvwEndData
            Tvwobj += obj_exp_datum(sol, obj; verbose=verbose)
        elseif obj isa EndTimeData
            tobj += obj_exp_datum(sol, obj; verbose=verbose)
        end
    end
    verbose && @info "loss call" tobj Tfobj Tvwobj tweight
    return ustrip(u"K^2", Tfobj + Tvw_weight*Tvwobj + tweight*tobj)
end
obj_exp(sol::Val{NaN}, efd; kwargs...) = Inf

"""
    $(SIGNATURES)

A thin wrapper on [`obj_exp`](@ref), for backwards compatibility.
"""
function obj_expT(sol, efd;
    tweight=1.0u"K^2/hr^2", verbose = false, Tvw_weight=1.0)
    return obj_exp(sol, efd; tweight, verbose, Tvw_weight)
end

# --- Per-data-type residual counts ----------------------------------------

num_errs(obj::TfData) = length(obj.t_range)
num_errs(obj::TvwSeriesData) = length(obj.t_range)
num_errs(obj::TvwEndData) = 1
num_errs(obj::EndTimeData) = 1

"""
    $(SIGNATURES)

Compute the number of data points available in `efd` for comparison to model solution.

This is useful for caching a residual vector for least-squares fitting, e.g. with [`err_exp!`](@ref) and [`nls_pd!`](@ref).
"""
function num_errs(efd::ExpFitData)
    nerr = mapreduce(num_errs, +, efd.data, init=0)
    return nerr
end

# --- Per-data-type residual functions ---------------------------------------
# Each fills `errs[i0+1 : i0+n]` with the residuals for one data object,
# returning the number of slots consumed. `not_avail_err` fills points where
# the model solution is unavailable (before `i_solstart` or after `trim`).

const not_avail_err = 0.0 # an error value to return for points where the solution is unavailable, e.g. if the model dries faster

function err_Tf!(errs, i0, sol::ODESolution, obj::TfData, st::SolTrim; verbose=false)
    Tmd = model_result(sol, st, 2; verbose=verbose)
    Tf = obj.Tf
    itf = length(obj.t_range)
    trim = min(itf, length(Tmd))
    Tferrs = (Tf[st.i_solstart:trim] .- Tmd[begin:trim-st.i_solstart+1])/sqrt(trim-st.i_solstart+1)
    errs[i0+1:i0+st.i_solstart] .= not_avail_err
    errs[i0+st.i_solstart:i0+trim] .= ustrip.(u"K", Tferrs)
    errs[i0+trim+1:i0+itf] .= not_avail_err
    return itf
end

function err_Tvw!(errs, i0, sol::ODESolution, obj::TvwSeriesData, st::SolTrim; verbose=false)
    Tvwmd = model_result(sol, st, 3; verbose=verbose)
    Tvw = obj.Tvw
    itvw = length(obj.t_range)
    trim = min(itvw, length(Tvwmd))
    Tvw_errs = (Tvw[st.i_solstart:trim] .- Tvwmd[begin:trim-st.i_solstart+1])/sqrt(trim-st.i_solstart+1)
    errs[i0+1:i0+st.i_solstart] .= not_avail_err
    errs[i0+st.i_solstart:i0+trim] .= ustrip.(u"K", Tvw_errs)
    errs[i0+trim+1:i0+itvw] .= not_avail_err
    return itvw
end

function err_Tvw!(errs, i0, sol::ODESolution, obj::TvwEndData; verbose=false)
    Tvw_err = sol[3, end]*u"K" - uconvert(u"K", obj.Tvw_end)
    errs[i0+1] = ustrip(u"K", Tvw_err)
    return 1
end

function err_tend!(errs, i0, sol::ODESolution, obj::EndTimeData; tweight=1.0u"K/hr", verbose=false)
    tmd = sol.t[end]*u"hr"
    t_end = obj.t_end
    if t_end isa Tuple # See if is inside window and scale appropriately
        mid_t = (t_end[1] + t_end[2]) / 2.0
        if tmd < t_end[1]
            t_err = (mid_t - tmd)
        elseif tmd > t_end[2]
            t_err = (mid_t - tmd)
        else # Inside window, so no error
            t_err = 0.0u"hr"
        end
    else
        t_err = (t_end - tmd)
    end
    errs[i0+1] = ustrip(u"K", t_err*tweight)
    return 1
end

const errexp_doc = """
Evaluate the error between model solution `sol` and experimental data in `efd`.

The in-place version `err_exp!` fills the `errs` vector with the errors, while the non-in-place version `err_exp` returns a new vector of errors.
In-place `err_exp!` thus doesn't allocate, but requires `errs` to have length `num_errs(efd)`.

In contrast to `obj_exp()`, which sums all the squared residuals, this function fills the passed array
`errs` with each separate residual, which is suited for least squares algorithms.
- `errs` is a vector of length `num_errs(efd)`, which this function fills with the errors.
-`sol` is a solution to an appropriate model; see [`gen_sol_pd`](@ref) for a helper function.
- `efd` is an instance of [`ExpFitData`](@ref), which contains some information about what to compare.
- `tweight = 1.0u"K/hr"` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
Each time series, plus the end time, is given equal weight by dividing by its length; error is given in K (but `ustrip`ped).

Note that if `efd` has vial wall temperatures (i.e. a [`TvwSeriesData`](@ref) or [`TvwEndData`](@ref) object in `efd.data`), the third-index variable in `sol` is assumed to be temperature, as is true for solutions with [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `efd`, error is computed for each separately; likewise for `Tvw`.
"""

"""
    $(SIGNATURES)

$errexp_doc
"""
function err_exp!(errs, sol::ODESolution, efd; tweight=1.0u"K/hr", verbose = false)
    if length(errs) != num_errs(efd)
        error("Wrong length of cached residual vector.")
    end
    # errs .= 0.0 # If indexing is handled correctly, this should not be necessary.
    if sol.retcode != ReturnCode.Terminated || length(sol.u) <= 2
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
        errs .= Inf
        return
    end
    # Check which, if any, of the data objects have their own time vector provided
    separate_time = map((x->has_timevec(x) && !ismissing(x.t)), efd.data)
    container_trim = trim_sol(sol, efd.t)
    last_ind = 0
    for (obj, sep_trim) in zip(efd.data, separate_time)
        st = if sep_trim
            trim_sol(sol, fit_t(efd, obj))
        elseif has_timevec(obj)
            container_trim
        else
            nothing
        end
        if obj isa TfData
            last_ind += err_Tf!(errs, last_ind, sol, obj, st; verbose)
        elseif obj isa TvwSeriesData
            last_ind += err_Tvw!(errs, last_ind, sol, obj, st; verbose)
        elseif obj isa TvwEndData
            last_ind += err_Tvw!(errs, last_ind, sol, obj; verbose)
        elseif obj isa EndTimeData
            last_ind += err_tend!(errs, last_ind, sol, obj; tweight, verbose)
        end
    end
    if last_ind != length(errs)
        error("Indexing problems...")
    end
    verbose && @info "loss call" sol.t[end]*u"hr" size(errs) sum(abs2.(errs))
    return nothing
end
err_exp!(errs, sol::Val{NaN}, efd; kwargs...) = errs .= Inf;

"""
    $(SIGNATURES)

$errexp_doc
"""
function err_exp(sol, efd; tweight=1.0u"K/hr", verbose = false)
    errs = zeros(num_errs(efd))
    err_exp!(errs, sol, efd; tweight, verbose)
    return errs
end
