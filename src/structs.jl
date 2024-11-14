export RpFormFit, RampedVariable, ConstPhysProp, PrimaryDryFit

"""
A convenience type for dealing with the common functional form given to Rp and Kv.

An object `Rp = RpFormFit(A, B, C)` can be called as `Rp(x)`, which simply computes `A + B*x/(1 + C*x)`.
Likewise, `Kv = RpFormFit(Kc, Kp, Kd)` can be called as `Kv(p)` to get `Kc + Kp*p/(1 + Kd*p)`.

Be careful to pass dimensionally consistent values.
"""
struct RpFormFit{T1, T2, T3}
    R0::T1
    A1::T2
    A2::T3
end

function (ff::RpFormFit)(x)
    return (ff.R0 + ff.A1*x/(1+ustrip(NoUnits, ff.A2*x)))
end

function Base.show(io::IO, rp::RpFormFit) 
    return print(io, "RpFormFit($(rp.R0), $(rp.A1), $(rp.A2))")
end

"""
A convenience type for computing temperatures, pressures, etc. with multiple setpoints in sequence,
and linear interpolation according to a fixed ramp rate between set points

Three main constructors are available:
For a non-varying value, call with one argument:

    RampedVariable(constant_setpt)

For one ramp from initial value to set point with indefinite hold, call with two arguments:

    RampedVariable(setpts, ramprate)
    
And for multiple setpoints, call with three arguments:

    RampedVariable(setpts, ramprates, holds)

With three arguments, `setpts`, `ramprates`, and `holds` should all be vectors, with lengths N+1, N, N-1 respectively.

The resulting RampedVariable `rv = RampedVariable(...)` can be called as `rv(x)` at any (dimensionally consistent) value of x, 
and will return the value at that time point along the ramp process.

A plot recipe is also provided for this type, e.g. `plot(rv; tmax=10u"hr")` where `tmax` indicates where to stop drawing the last setpoint hold.
"""
struct RampedVariable{vary, T1, T2, T3, T4}
    setpts::Union{T1, Vector{T1}}
    ramprates::Union{Vector{T2}, T2}
    holds::Union{Vector{T3}, T3}
    timestops::Union{Vector{T4}, T4}
end

# get_dimensions(::Type{Quantity{T,D,U}}) where {T, D, U} = D
function (rv::RampedVariable{false, T1,T2,T3,T4})(t) where {T1,T2,T3,T4}
    return rv.setpts::T1
end
function (rv::RampedVariable{true, T1,T2,T3,T4})(t) where {T1,T2,T3,T4}
    im = findlast(rv.timestops .<= t)
    if im == length(rv.timestops)
        return rv.setpts[end]
    elseif iseven(im)
        return rv.setpts[im÷2+1]
    else
        ip = im+1
        return ((rv.setpts[ip÷2+1] - rv.setpts[ip÷2])/(rv.timestops[ip] - rv.timestops[im])*(t - rv.timestops[im]) + rv.setpts[ip÷2])
    end
end


function RampedVariable(hold)
    RampedVariable{false, typeof(hold), Nothing, Nothing, Nothing}(hold, nothing, nothing, nothing)
end

function RampedVariable(setpts, ramprate)

    if length(ramprate) == 0 || length(setpts) == 1
        @error "If no ramp necessary, construct RampedVariable with only one argument." ramprate
    end
    if length(ramprate) >= 2 || length(setpts) > 2
        @error "For multiple ramps, need at least one hold time. Construct RampedVariable with three arguments." ramprate
    end
    if length(setpts) != 2
        @error "Number of set points should be 1 more than ramps, since initial is included"
    end
    timestops = fill(0.0*setpts[1]/ramprate[1], 2)
    timestops[2] = timestops[1] + (setpts[2]-setpts[1])/ramprate
    RampedVariable{true, eltype(setpts), typeof(ramprate), Nothing, eltype(timestops)}(setpts, [ramprate], nothing, timestops)
end

function RampedVariable(setpts, ramprates, holds)

    if (length(ramprates) != length(holds) + 1 ) || (length(ramprates)==0)
        @error "Number of ramps should be zero or number of holds + 1"
    end
    if length(setpts) != length(ramprates) + 1
        @error "Number of set points should be 1 more than ramps, since initial is included"
    end
    timestops = fill(0.0*setpts[1]/ramprates[1], length(ramprates) + length(holds) + 1)
    (ramp, rest) = Iterators.peel(ramprates)
    timestops[2] = timestops[1] + (setpts[2]-setpts[1])/ramp
    if timestops[2] < timestops[1]
        @warn "Ramp rate given with probably the wrong sign, changing its sign"
        timestops[2] = timestops[1] + (timestops[1]-timestops[2])
    end
    for (i, ramp) in enumerate(rest)
        timestops[2i+1] = timestops[2i] + holds[i]
        timestops[2i+2] = timestops[2i+1] + (setpts[i+2]-setpts[i+1])/ramp
        if timestops[2i+2] < timestops[2i+1]
            @warn "Ramp rate given with probably the wrong sign, changing its sign"
            timestops[2i+2] = timestops[2i+1] + (timestops[2i+1]-timestops[2i+2])
        end
    end
    RampedVariable{true, eltype(setpts), eltype(ramprates), eltype(holds), eltype(timestops)}(setpts, ramprates, holds, timestops)
end

function Base.hash(rv::RampedVariable, h::UInt)
    hash(rv.setpts, hash(rv.ramprates, hash(rv.holds, hash(rv.timestops, hash(:RampedVariable, h)))))
end

function Base.show(io::IO, rv::RampedVariable{false, T1,T2,T3,T4}) where {T1,T2,T3,T4}
    return print(io, "RampedVariable($(rv.setpts))")
end
    
function Base.show(io::IO, rv::RampedVariable{true, T1,T2,T3,T4}) where {T1,T2,T3,T4}
    if length(rv.setpts) == 2
        return print(io, "RampedVariable($(rv.setpts), $(rv.ramprates[1]))")
    else 
        return print(io, "RampedVariable($(rv.setpts), $(rv.ramprates), $(rv.holds))")
    end
end

struct PhysProp{V, Tdep, Pdep, Fdep}
    valfunc::V
end

(pp::PhysProp{V, true , true , true } where V)(T, p, f) = pp.valfunc(T, p, f)
(pp::PhysProp{V, false, true , true } where V)(T, p, f) = pp.valfunc(p, f)
(pp::PhysProp{V, true , false, true } where V)(T, p, f) = pp.valfunc(T, f)
(pp::PhysProp{V, true , true , false} where V)(T, p, f) = pp.valfunc(T, p)
(pp::PhysProp{V, false, false, true } where V)(T, p, f) = pp.valfunc(f)
(pp::PhysProp{V, true , false, false} where V)(T, p, f) = pp.valfunc(T)
(pp::PhysProp{V, false, true , false} where V)(T, p, f) = pp.valfunc(p)
(pp::PhysProp{V, false, false, false} where V)(T, p, f) = pp.valfunc

PhysProp(x) = PhysProp{typeof(x), false, false, false}(x::T)
function PhysProp(x, args...) 
    func_T = (:T ∈ args)
    func_p = (:p ∈ args)
    func_f = (:f ∈ args)
    PhysProp{typeof(x), func_T, func_p, func_f}(x)
end

struct ConstPhysProp 
    val
end
(cpp::ConstPhysProp)(args...) = cpp.val

# function Base.show(io::IO, pp::PhysProp) 
#     return print(io, "PhysProp($(rv.setpts[1]))")
# end

"""
PrimaryDryFit: a type for storing experimental data and indicating how it should be fit.

Provided constructors:

    PrimaryDryFit(t_Tf, Tfs, Tf_iend, t_Tvw, Tvws, Tvw_iend, t_end)
    PrimaryDryFit(t_Tf, Tf::V) where V<:AbstractVector
    PrimaryDryFit(t_Tf, Tfs::T) where T<:Tuple
    PrimaryDryFit(t_Tf, Tfs, t_end::T) where T<:Unitful.Time
    PrimaryDryFit(t_Tf, Tfs, Tvw_end::T) where T<:Unitful.Temperature
    PrimaryDryFit(t_Tf, Tfs, Tvw_end::T1, t_end::T2) where {T1<:Unitful.Temperature, T2<:Unitful.Time}
    PrimaryDryFit(t_Tf, Tfs, t_Tvw, Tvws) 
    PrimaryDryFit(t_Tf, Tfs, t_Tvw, Tvws, t_end)

The use of this struct is determined in large part by the implementation of 
[`LyoPronto.obj_expcomp`](@ref). If a given field is not available, set it
to `missing` and things should basically work. At least `t_Tf` and `Tfs` are 
expected to always be provided.

With the exception of the two-argument `(t_Tf, Tf)` constructor, `Tfs` and `Tvws` should always be tuples of vectors (one vector per time series).

`Tf_iend` and `Tvw_iend` default to `[length(Tf) for Tf in Tfs]` and `[length(Tvw) for Tvw in Tvws]`, 
respectively, with one value for each temperature series; 
they are used to dictate if a given temperature series should
be truncated sooner than the full length in the fitting procedure.
This implies that all the `Tf` temperature series are valid initially at the same
time points, then stop having measured values after a different number of measurements.

`t_end` indicates an end of drying, particularly if taken from other measurements
(e.g. from Pirani-CM convergence). If set to `missing`, it is ignored in the
objective function.

Principal Cases:
- Conventional: provide only `t_Tf, Tfs`
- RF with measured vial wall: provide `t_Tf, Tfs, t_Tvw, Tvws`, 
- RF, matching model Tvw to experimental Tf[end] without measured vial wall: provide `t_Tf, Tfs, Tvws`, set `t_Tvw` to `missing`
"""
struct PrimaryDryFit_7{T1, T2, T3, T4} 
    t_Tf::Vector{T1}
    Tfs::Tuple{Vector{T2}, Vararg{Vector{T2}}}
    Tf_iend::Vector{T3}# = [length(Tf) for Tf in Tfs]
    t_Tvw::Union{Missing, Vector{T1}}# = missing
    Tvws::Union{Missing, T2, Tuple{Vector{T2}, Vararg{Vector{T2}}}}# = missing
    Tvw_iend::Union{Missing, Vector{T3}}# = (ismissing(Tvws) ? missing : [length(Tvw) for Tvw in Tvws])
    t_end::T4# = missing
end

# Primary constructor
function PrimaryDryFit_7(t_Tf, Tfs, t_Tvw, Tvws, t_end) 
    if ismissing(t_Tvw) && !ismissing(Tvws) && (length(Tvws) > 1)
        throw("If no time passed for `t_Tvw`, `Tvws` should be a single value, treated as a vial endpoint temperature")
    end
    PrimaryDryFit_7(t_Tf, Tfs, [length(Tf) for Tf in Tfs], t_Tvw, Tvws, 
    ((ismissing(Tvws) || ismissing(t_Tvw)) ? missing : [length(Tvw) for Tvw in Tvws]), t_end)
end
# Convenience constructors
PrimaryDryFit_7(t_Tf, Tfs, t_Tvw, Tvws)  = PrimaryDryFit_7(t_Tf, Tfs, t_Tvw, Tvws, missing)
PrimaryDryFit_7(t_Tf, Tfs, Tvw::Unitful.Temperature)  = PrimaryDryFit_7(t_Tf, Tfs, missing, Tvw, missing)
PrimaryDryFit_7(t_Tf, Tfs, Tvw::Unitful.Temperature, t_end::Unitful.Time)  = PrimaryDryFit_7(t_Tf, Tfs, missing, Tvw, t_end)
PrimaryDryFit_7(t_Tf, Tfs, t_end::Unitful.Time)  = PrimaryDryFit_7(t_Tf, Tfs, missing, missing, t_end)
PrimaryDryFit_7(t_Tf, Tf::V) where V<:Vector = PrimaryDryFit_7(t_Tf, (Tf,), missing, missing, missing)
PrimaryDryFit_7(t_Tf, Tfs::T) where T<:Tuple = PrimaryDryFit_7(t_Tf, Tfs, missing, missing, missing)

# Temporary kludge to make Revise work better
PrimaryDryFit = PrimaryDryFit_7

function Base.:(==)(p1::PrimaryDryFit, p2::PrimaryDryFit)
    cond1 = p1.t_Tf == p2.t_Tf
    cond2 = p1.Tfs == p2.Tfs
    cond3 = p1.Tf_iend == p2.Tf_iend
    cond4 = ismissing(p1.t_Tvw) ? ismissing(p2.t_Tvw) : (p1.t_Tvw == p2.t_Tvw)
    cond5 = ismissing(p1.Tvws) ? ismissing(p2.Tvws) : (p1.Tvws == p2.Tvws)
    cond6 = ismissing(p1.Tvw_iend) ? ismissing(p2.Tvw_iend) : (p1.Tvw_iend == p2.Tvw_iend)
    cond7 = ismissing(p1.t_end) ? ismissing(p2.t_end) : (p1.t_end == p2.t_end)
    return all([cond1, cond2, cond3, cond4, cond5, cond6, cond7])
end