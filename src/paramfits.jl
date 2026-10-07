
# Data objects for ExpFitData. Each holds one kind of experimental measurement
# plus an optional per-object time vector (`t`). When `t` is `missing`, the
# container's `t` is used (see `fit_t`).
"""
    AbstractExpDatum

An abstract type for all types of experimental data used in fitting.

If you add a new type with this abstract type, you should define the following methods for your type:
- `resid_name(::MyDatum)::Symbol`: return a unique symbol identifying the type of residuals associated with your data type. This is used to map to a weight in the `weights` NamedTuple passed to [`obj_exp`](@ref) and [`err_exp`](@ref).
- `time_bound_data(::MyDatum)`: return `true` if your data type is associated with a time vector, and `false` otherwise. If `true`, you should also define:
  - `has_timevec(md::MyDatum)`: return `true` if the instance `md` has its own time vector, and `false` if it uses the container's time vector.
  - `nontrivial_t_range(md::MyDatum)`: return `true` if the instance `md` has a nontrivial set of time indices (i.e., not starting at 1 and incrementing by 1) into the container's or instance's time vector, and `false` otherwise.

"""
abstract type AbstractExpDatum end

"""
    time_bound_data(obj::AbstractExpDatum)::Bool
Trait for identifying whether a given type should be associated with a time vector
"""
function time_bound_data end

"""
    $(SIGNATURES)
Check if a given instance of an `AbstractExpDatum` subtype has its own self-contained time vector.

This is strictly false if `time_bound_data` is false
"""
function has_timevec end

"""
    $(SIGNATURES)
Check if a given instance of an `AbstractExpDatum` subtype has a nontrivial set of time indices.
    
"Nontrivial" in this case means that the indices do not start at 1 and/or do not increment by 1, i.e. the data series is not aligned with the start of the time vector or is not contiguous.
"""
function nontrivial_t_range end

"""
resid_name(obj::AbstractExpDatum)::Symbol

Return a symbol identifying the type of residuals associated with `obj`. This is used to map to a weight in the `weights` NamedTuple passed to [`obj_exp`](@ref) and [`err_exp`](@ref). The default names are:
- `:Tf` for [`TfData`](@ref)
- `:Tvw` for [`TvwSeriesData`](@ref) and [`TvwEndData`](@ref)
- `:endTime` for [`EndTimeData`](@ref)

If you define your own type of experimental data, you must define `resid_name(::MyDatum)::Symbol`
to return a unique symbol for your data type, and then map that symbol to a weight in the 
`weights` NamedTuple passed to [`obj_exp`](@ref) and [`err_exp`](@ref).
"""
function resid_name end

"""
    TfData(Tf; t=missing)
    TfData(Tf, t_range; t=missing)

A single experimental frozen-product (Tf) temperature series, for use in [`ExpFitData`](@ref).

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

Base.length(td::TfData) = length(td.t_range)
time_bound_data(::TfData) = true
has_timevec(td::TfData) = !ismissing(td.t)
nontrivial_t_range(td::TfData) = first(td.t_range) > 1 || step(td.t_range) != 1
resid_name(::TfData) = :Tf

# Provide separate types for Tvw as a full series or as a single endpoint
"""
    TvwSeriesData(Tvws; t=missing)
    TvwSeriesData(Tvws, t_range; t=missing)

A single experimental vial-wall (Tvw) temperature series, for use in [`ExpFitData`](@ref).

Fields:
- `Tvw`: an `AbstractVector` of `Temperature` (one series).
- `t_range`: a range of integer indices into the time vector that `Tvw` corresponds to;
  `Tvw[k]` is the temperature measured at time `t[t_range[k]]`.
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

Base.length(td::TvwSeriesData) = length(td.t_range)
time_bound_data(::TvwSeriesData) = true
has_timevec(td::TvwSeriesData) = !ismissing(td.t)
nontrivial_t_range(td::TvwSeriesData) = first(td.t_range) > 1 || step(td.t_range) != 1
resid_name(::TvwSeriesData) = :Tvw

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

time_bound_data(::TvwEndData) = false
has_timevec(::TvwEndData) = false
resid_name(::TvwEndData) = :Tvw

"""
    EndTimeData(t_end)

An experimentally measured end of primary drying, for use in [`ExpFitData`](@ref).

`t_end` is typically determined from non-temperature measurements (e.g. Pirani-CM
convergence). It can be a single time, or a tuple of two times indicating a window of
acceptable drying times: in the objective function, any model drying time inside that
window is not penalized, while outside the window a squared error applies, as for the
single-time case. If no [`EndTimeData`](@ref) object is present in the fit, the end of
drying is ignored in the objective function.
"""
struct EndTimeData{T1} <: AbstractExpDatum
    t_end::T1      # Unitful.Time OR Tuple of two Times
    function EndTimeData(t_end)
        if !(t_end isa Unitful.Time) && !(t_end isa Tuple && length(t_end) == 2 && all(x -> x isa Unitful.Time, t_end))
            throw(ArgumentError("t_end should be a time or a tuple of two times"))
        end
        if t_end isa Tuple
            t_end = extrema(t_end) # Reorder to be (min, max)
        end
        new{typeof(t_end)}(t_end)
    end
end

time_bound_data(::EndTimeData) = false
has_timevec(::EndTimeData) = false
resid_name(::EndTimeData) = :t

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

Return the time vector to use for `obj`.
    
Three cases:
- `obj.t[obj.t_range]` if `obj.t` is set 
- `fitdat.t[obj.t_range]` if `t_range` makes sense and `obj.t` is `missing`
- `fitdat.t` if `obj` is not tied to a time vector

"""
function fit_t(fitdat::ExpFitData, dat::AbstractExpDatum)
    if time_bound_data(dat) 
        return has_timevec(dat) ? dat.t[dat.t_range] : fitdat.t[dat.t_range] 
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

Each data object is documented with its own fields and constructors. In particular, the
temperature series ([`TfData`](@ref), [`TvwSeriesData`](@ref)) carry a `t_range` selecting the
time points they correspond to, and any data object may carry its own time vector `t` (see
[`fit_t`](@ref)). An end of drying is supplied via an [`EndTimeData`](@ref) object and is
ignored in the objective function if absent.

Currently-implemented data objects include:
- [`TfData`](@ref): a single frozen-product temperature series
- [`TvwSeriesData`](@ref): a single vial-wall temperature series
- [`TvwEndData`](@ref): a single vial-wall temperature endpoint
- [`EndTimeData`](@ref): a single end-of-drying time, or a tuple of two times indicating an acceptable window

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
        tvw_obj = Tvws isa Number ? Tuple(TvwEndData(Tvws)) : Tuple(TvwSeriesData(Tvw) for Tvw in Tvws)
        objs = (objs..., tvw_obj...)
    end
    if !ismissing(t_end)
        objs = (objs..., EndTimeData(t_end))
    end
    data = objs
    ExpFitData(t, data)  # calls the actual current constructor
end
PrimaryDryFit(t, Tfs, Tvws, t_end) = PrimaryDryFit(t, Tfs; Tvws=Tvws, t_end=t_end)


function Base.:(==)(a::TfData, b::TfData)
    return a.Tf == b.Tf && a.t_range == b.t_range && isequal(a.t, b.t)
end
function Base.:(==)(a::TvwSeriesData, b::TvwSeriesData)
    return a.Tvw == b.Tvw && a.t_range == b.t_range && isequal(a.t, b.t)
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

"""
    SolTrim

A struct to contain information about which model points overlap with experiment.

If `preinterp` is true, then the model solution is pre-interpolated to the experimental time points;
the model's time indices `1:(length(sol.t)-1)` correspond exactly to the experimental
time indices `ti_model_start:ti_model_end`. 

If `preinterp` is false, then the model solution is not pre-interpolated, and the model solution will need to be interpolated to the experimental time points.
at time indices `ti_m_start:ti_m_end`. 
"""
struct SolTrim{I <: Integer, T}
    preinterp::Bool
    ti_m_start::I
    ti_m_end::I
    len::I
    t::T
end

"""
    $(SIGNATURES)
Return the range of time indices for the experimental time points that match up with the model solution.
"""
function exp_time_inds(st::SolTrim)
    return st.ti_m_start:st.ti_m_end
end
"""
    $(SIGNATURES)
Return the range of time indices for the model solution that match up with the experimental time points.
"""
function model_time_inds(st::SolTrim)
    return 1:st.len
end
Base.length(st::SolTrim) = st.len


function check_preinterp(sol::ODESolution, t_exp)
    nt = min(length(sol.t) - 1, length(t_exp)) # possible number of valid times
    ti_m_start = searchsortedfirst(t_exp, sol.t[begin]*u"hr")
    # Check if the solution is pre-interpolated to the time points in t
    preinterp = mapreduce(≈, &, sol.t[1:nt]*u"hr", t_exp[ti_m_start:ti_m_start+nt-1])
    return preinterp, ti_m_start
end

"""
    trim_sol(sol::ODESolution, t_exp)

Trim the solution `sol` to the experimental time grid `t_exp`, returning a `SolTrim` object.

As a first pass, this should be called on an `ExpFitData`.
Thenshould be called separately for data objects which have their own time vector.
"""
function trim_sol(sol::ODESolution, t_exp, preinterp::Bool, ti_m_start)
    # Last index of t_exp that is <= sol.t[end]
    nt = min(length(sol.t) - 1, length(t_exp))
    ti_m_end = searchsortedlast(t_exp, sol.t[end]*u"hr")

    # Make sure there is at least one time point overlapping
    if ti_m_end - ti_m_start < 1
        throw(ArgumentError("Experimental time vector does not overlap with model solution time points"))
    end
        
    #TODO: see if we can avoid allocating a new time vector here, for performance
    # First idea of checking of passing `missing` for preinterpolated is type-unstable, I think
    # pass_t = preinterp ? t_exp[1:0] : t_exp[ti_m_start:ti_m_end]
    return SolTrim(preinterp, ti_m_start, ti_m_end, ti_m_end - ti_m_start + 1, t_exp[ti_m_start:ti_m_end])
end

"""
    $(SIGNATURES)

Compute the model result for a given solution `sol` with time trimming `st`, variable `idx`.

If `idx` is specified, the function evaluates the solution at that index and gives it the 
units specified by the `unit` kwarg (defaulting to `u"K"`). 

If `idx` is not specified, the function evaluates [`calc_md_Q`](@ref) at all appropriate 
time points and returns a corresponding Table. If only a single column from that table is desired,
pass the keyword argument `var` to select the column by name (e.g. `var=:md`).

If you have a choice between the two, it will be more efficient to use `idx` since it only
has to interpolate the solution for that index, while the `var` option will evaluate the full model at all time points and then select the column.
"""
function model_result(sol::ODESolution, st::SolTrim, idx; unit=u"K", verbose=false)
    if st.preinterp
        # We want multi-D array indexing, so index into sol, not sol.u
        return sol[idx, begin:end-1]*unit # Leave off last time point because is end time
    end
    res = sol.(ustrip.(u"hr", st.t), idxs=idx)*unit
    # Sometimes the interpolation procedure of the solution produces temperatures below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if first(res) isa Unitful.Temperature && any(res .< 0u"K")
        subzero = findall(Vector(res .< 0u"K"))
        res[subzero] .= res[subzero[1] - 1]
        verbose && @info "bad interpolation" subzero res[subzero]
    end
    return res
end
function model_result(sol::ODESolution, st::SolTrim; var=nothing, verbose=false)
    t_trim_nd = ustrip.(u"hr", st.t)
    uu = if st.preinterp
        # We want a vector of vectors, so access sol.u
        sol.u[begin:end-1] # Leave off last time point because is end time
    else
        uu = sol.(t_trim_nd)
    end
    if isnothing(var)
        mdq = map((u, t) -> calc_md_Q(u, sol.prob.p, t), uu, t_trim_nd)
        return Table(mdq)
    else
        varq = map((u, t) -> calc_md_Q(u, sol.prob.p, t)[var], uu, t_trim_nd)
        return varq
    end
end
# TODO: decide whether to provide the following. Probably not worth it, since in the absence
# of a SolTrim the indexing is pretty trivial.
# function model_result(sol::ODESolution, idx; unit=u"K", verbose=false)
#     res = sol[idx, :]*unit 
#     if first(res) isa Unitful.Temperature && any(res .< 0u"K")
#         subzero = findall(Vector(res .< 0u"K"))
#         res[subzero] .= res[subzero[1] - 1]
#         verbose && @info "bad interpolation" subzero res[subzero]
#     end
#     return res
# end
# function model_result(sol::ODESolution; var=nothing, verbose=false)
#     uu = sol.u     
#     if isnothing(var)
#         mdq = map((u, t) -> calc_md_Q(u, sol.prob.p, t), uu, sol.t)
#         return Table(mdq)
#     else
#         varq = map((u, t) -> calc_md_Q(u, sol.prob.p, t)[var], uu, t_trim_nd)
#         return varq
#     end
# end

# ---------------
# Weighting functions for squared-error loss and for residuals. These are used in `obj_exp` and `err_exp`, respectively.

"""
    $(SIGNATURES)

Construct a NamedTuple of weights per experiment type for use with [`obj_exp`](@ref).
    
`weights` maps the names returned by [`resid_name`](@ref) to inverse-squared-unit weights, 
so even custom experimental data types can be weighted against each other without changing 
LyoPronto. `t`, `Tf`, and `Tvw` have built-in defaults.
    
The names of each kwarg must match the results of [`resid_name`](@ref) for each type of 
`AbstractExpDatum`, e.g. `:Tf`, `:t`, `:Tvw` for the builtin types.

The value of each kwarg should have inverse-square units, e.g. `u"hr^-2"` for `:t`.
"""
function loss_weighting(;
    t = 1.0u"hr^-2", 
    Tvw = 1.0u"K^-2", 
    Tf = 1.0u"K^-2",
    kwargs...
    )
    return (;t, Tvw, Tf, kwargs...)
end

"""
    $(SIGNATURES)

Construct a NamedTuple of weights per experiment type for use with [`err_exp`](@ref).
    
`weights` maps the names returned by [`resid_name`](@ref) to inverse-unit weights, 
so even custom experimental data types can be weighted against each other without changing 
LyoPronto. `t`, `Tf`, and `Tvw` have built-in defaults.
    
The names of each kwarg must match the results of [`resid_name`](@ref) for each type of 
`AbstractExpDatum`, e.g. `:Tf`, `:t`, `:Tvw` for the builtin types.

The value of each kwarg should have inverse units, e.g. `u"hr^-1"` for `:t`.
"""
function residual_weighting(;
    t = 1.0u"hr^-1", 
    Tvw = 1.0u"K^-1", 
    Tf = 1.0u"K^-1",
    kwargs...
    )
    return (;t, Tvw, Tf, kwargs...)
end

# --- Per-data-type objective functions ------------------------------------
# Each returns a scalar contribution (in K^2 for temperature, hr^2 for time).

function obj_Tf(sol::ODESolution, obj::TfData, t; verbose=false)
    st = trim_sol(sol, t)
    obj_Tf(sol, st, obj; verbose)
end
function obj_Tf(sol::ODESolution, st::SolTrim, dat::TfData; verbose=false)
    Tmd = model_result(sol, st, 2; verbose) # Tf at index 2
    verbose && @info "Tf_model = $Tmd"
    resid = sum(abs2, (dat.Tf[exp_time_inds(st)] .- Tmd[model_time_inds(st)]))/(length(st))
    verbose && @info "Tf_err = $resid"
    return resid
end

function obj_Tvw(sol::ODESolution, obj::TvwSeriesData, t; verbose=false)
    st = trim_sol(sol, t)
    obj_Tvw(sol, st, obj; verbose)
end
function obj_Tvw(sol::ODESolution, st::SolTrim, dat::TvwSeriesData; verbose=false)
    Tmd = model_result(sol, st, 3; verbose) # Tvw at index 3
    resid = sum(abs2, (dat.Tvw[exp_time_inds(st)] .- Tmd[model_time_inds(st)]))/(length(st))
    verbose && @info "Tvw_err = $resid"
    return resid
end

function obj_Tvw(sol::ODESolution, obj::TvwEndData; verbose=false)
    resid = (sol[3, end]*u"K" - uconvert(u"K", obj.Tvw_end))^2
    verbose && @info "Tvw_err = $resid"
    return resid
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
- `weights` is a `NamedTuple` maps each experimental datum to an inverse-squared-unit weight;
    the default is constructed with [`loss_weighting`](@ref)().

To add a custom type of experimental data, 
define `resid_name(::MyDatum)::Symbol` and map that symbol to a weight accessed in a `NamedTuple`,
e.g. `weights[:t]` should be e.g. `1.0u"hr^-2"`.

Note that if `efd` has vial wall temperatures (i.e. a [`TvwSeriesData`](@ref) or [`TvwEndData`](@ref) object in `efd.data`), the third-index variable in `sol` is assumed to be temperature, as is true for the lumped capacitance model (see [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `efd`, squared error is computed for each separately then summed; likewise for `Tvw`.
"""
function obj_exp(sol::ODESolution, efd::ExpFitData;
    weights = loss_weighting(), verbose = false)
    if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
        verbose && @warn "ODE solve did not reach end of drying. Either parameters are bad, or tspan is not large enough." sol.retcode sol.prob.p.hf0 sol[end]
        return Inf
    end
    # Check if the 
    preinterp_shared, ti_m_start_shared = check_preinterp(sol, efd.t)
    verbose && @info "loss call" # Verbose gets passed to separate calls
    obj = mapreduce(+, efd.data) do dat
        single_err = if time_bound_data(dat) 
            sep_trim = has_timevec(dat) || nontrivial_t_range(dat)
            ti_m_start = if sep_trim
                searchsortedfirst(fit_t(efd, dat), sol.t[begin]*u"hr")
            else
                ti_m_start_shared
            end
            # If the container t is preinterpolated, then specific is preintepolated exactly if it shares t
            # If the container t is not preinterpolated, the solution almost certainly isn't
            #   and it probably isn't worth checking
            preinterp = preinterp_shared ? ~sep_trim : false
            st = trim_sol(sol, fit_t(efd, dat), preinterp, ti_m_start)
            obj_exp_datum(sol, st, dat; verbose)
        else
            # If the data is not time-bound, then the solution doesn't need any trimming
            obj_exp_datum(sol, dat; verbose)
        end
        # Evaluate the residual, divide by weight matched to data type, 
        # check for nondimensional, then strip units
        return ustrip(NoUnits, single_err * weights[resid_name(dat)])
    end
    return obj
end
obj_exp(sol::Val{NaN}, efd; kwargs...) = Inf



"""
    $(SIGNATURES)

A thin wrapper on [`obj_exp`](@ref), for backwards compatibility.
"""
function obj_expT(sol, efd;
    tweight=1.0, verbose = false, Tvw_weight=1.0)
    return obj_exp(sol, efd; verbose,
    weights=loss_weighting(t=tweight*u"hr^-2", Tvw=Tvw_weight*u"K^-2"))
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

function err_Tf!(errs, i0, sol::ODESolution, dat::TfData, st::SolTrim, weight; verbose=false)
    Tmd = model_result(sol, st, 2; verbose=verbose)
    Tferrs = (dat.Tf[exp_time_inds(st)] .- Tmd[model_time_inds(st)])/sqrt(length(st))
    ntf = length(dat.Tf)
    # Broadcast across 1:ntf, not across exp_time_inds
    not_avail_inds = findall(1:ntf .∉ (exp_time_inds(st),))
    errs[i0 .+ not_avail_inds] .= not_avail_err
    errs[i0 .+ exp_time_inds(st)] .= ustrip.(NoUnits, Tferrs * weight)
    return ntf
end

function err_Tvw_series!(errs, i0, sol::ODESolution, dat::TvwSeriesData, st::SolTrim, weight; verbose=false)
    Tvwmd = model_result(sol, st, 3; verbose=verbose)
    Tvw_errs = (dat.Tvw[exp_time_inds(st)] .- Tvwmd[model_time_inds(st)])/sqrt(length(st))
    ntvw = length(dat.Tvw)
    # Broadcast across 1:ntvw, not across exp_time_inds
    not_avail_inds = findall(1:ntvw .∉ (exp_time_inds(st),)) # C
    ## Fill those with sentinel value, and fill the rest with the actual error
    errs[i0 .+ not_avail_inds] .= not_avail_err
    errs[i0 .+ exp_time_inds(st)] .= ustrip.(NoUnits, Tvw_errs * weight)
    return ntvw
end

function err_Tvw_end!(errs, i0, sol::ODESolution, obj::TvwEndData, weight; verbose=false)
    Tvw_err = sol[3, end]*u"K" - uconvert(u"K", obj.Tvw_end)
    errs[i0+1] = ustrip(NoUnits, Tvw_err * weight)
    return 1
end

function err_tend!(errs, i0, sol::ODESolution, obj::EndTimeData, weight; verbose=false)
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
    errs[i0+1] = ustrip(NoUnits, t_err*weight)
    return 1
end

function err_exp_datum!(errs, i0, sol::ODESolution, st::SolTrim, obj::TfData, weight; verbose=false)
    return err_Tf!(errs, i0, sol, obj, st, weight; verbose)
end
function err_exp_datum!(errs, i0, sol::ODESolution, obj::TfData, t, weight; verbose=false)
    st = trim_sol(sol, t)
    return err_Tf!(errs, i0, sol, obj, st, weight; verbose)
end
function err_exp_datum!(errs, i0, sol::ODESolution, st::SolTrim, obj::TvwSeriesData, weight; verbose=false)
    return err_Tvw_series!(errs, i0, sol, obj, st, weight; verbose)
end
function err_exp_datum!(errs, i0, sol::ODESolution, obj::TvwSeriesData, t, weight; verbose=false)
    st = trim_sol(sol, t)
    return err_Tvw_series!(errs, i0, sol, obj, st, weight; verbose)
end
function err_exp_datum!(errs, i0, sol::ODESolution, obj::TvwEndData, weight; verbose=false)
    return err_Tvw_end!(errs, i0, sol, obj, weight; verbose)
end
function err_exp_datum!(errs, i0, sol::ODESolution, obj::EndTimeData, weight; verbose=false)
    return err_tend!(errs, i0, sol, obj, weight; verbose)
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
- `weights` is a `NamedTuple` mapping each experimental datum to an inverse-unit weight; the default is constructed with [`residual_weighting`](@ref)`()`.
Each time series, plus the end time, is given equal weight by dividing by its length; error is given in K (but `ustrip`ped).

Note that if `efd` has vial wall temperatures (i.e. a [`TvwSeriesData`](@ref) or [`TvwEndData`](@ref) object in `efd.data`), the third-index variable in `sol` is assumed to be temperature, as is true for solutions with [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `efd`, error is computed for each separately; likewise for `Tvw`.
"""

"""
    $(SIGNATURES)

$errexp_doc
"""
function err_exp!(errs, sol::ODESolution, efd; weights=residual_weighting(), verbose = false)
    if length(errs) != num_errs(efd)
        error("Wrong length of cached residual vector.")
    end
    # errs .= 0.0 # If indexing is handled correctly, this should not be necessary.
    if sol.retcode != ReturnCode.Terminated || length(sol.u) <= 2
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
        errs .= Inf
        return
    end
    # Check which, if any, of the data objects have their own time vector or t_range provided
    preinterp_shared, ti_m_start_shared = check_preinterp(sol, efd.t)
    last_ind = 0
    for dat in efd.data
        last_ind += if time_bound_data(dat)
            sep_trim = has_timevec(dat) || nontrivial_t_range(dat)
            ti_m_start = if sep_trim
                searchsortedfirst(dat.t, sol.t[begin]*u"hr")
            else
                ti_m_start_shared
            end
            # If the container t is preinterpolated, then specific is preintepolated exactly if it shares t
            # If the container t is not preinterpolated, the solution almost certainly isn't
            #   and it probably isn't worth checking
            preinterp = preinterp_shared ? ~sep_trim : false
            st = trim_sol(sol, fit_t(efd, dat), preinterp, ti_m_start)
            err_exp_datum!(errs, last_ind, sol, st, dat, weights[resid_name(dat)]; verbose)
        else
            err_exp_datum!(errs, last_ind, sol, dat, weights[resid_name(dat)]; verbose)
        end
    end
    if last_ind != length(errs)
        error("Indexing problems. Should have reached the end of the residual vector, but last_ind = $last_ind, length(errs) = $(length(errs))")
    end
    verbose && @info "loss call" sol.t[end]*u"hr" size(errs) sum(abs2.(errs))
    return nothing
end
err_exp!(errs, sol::Val{NaN}, efd; kwargs...) = errs .= Inf;

"""
    $(SIGNATURES)

$errexp_doc
"""
function err_exp(sol, efd; weights=residual_weighting(), verbose = false)
    errs = zeros(num_errs(efd))
    err_exp!(errs, sol, efd; weights, verbose)
    return errs
end
