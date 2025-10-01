
@concrete terse struct RpFormFit
    R0
    A1
    A2
end
@doc """
A convenience type for dealing with the common functional form given to Rp and Kv.

An object `Rp = RpFormFit(A, B, C)` can be called as `Rp(x)`, which simply computes `A + B*x/(1 + C*x)`.
Likewise, `Kv = RpFormFit(Kc, Kp, Kd)` can be called as `Kv(p)` to get `Kc + Kp*p/(1 + Kd*p)`.

Be careful to pass dimensionally consistent values.
"""
RpFormFit

function (ff::RpFormFit)(x)
    return (ff.R0 + ff.A1*x/(1+ustrip(NoUnits, ff.A2*x)))
end

@concrete terse struct RampedVariable{vary}
    setpts
    ramprates
    holds
    timestops
end

@doc """
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
RampedVariable

# get_dimensions(::Type{Quantity{T,D,U}}) where {T, D, U} = D
function (rv::RampedVariable{false})(t)
    return rv.setpts
end
function (rv::RampedVariable{true})(t)
    im = searchsortedfirst(rv.timestops, t) - 1
    if im == 0 # Negative time
        return rv.setpts[1]
    elseif im == length(rv.timestops)
        return rv.setpts[end]
    elseif iseven(im)
        return rv.setpts[im÷2+1]
    else
        ip = im+1
        return ((rv.setpts[ip÷2+1] - rv.setpts[ip÷2])/(rv.timestops[ip] - rv.timestops[im])*(t - rv.timestops[im]) + rv.setpts[ip÷2])
    end
end


function RampedVariable(hold)
    RampedVariable{false}(hold, nothing, nothing, nothing)
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
    RampedVariable{true}(setpts, [ramprate], nothing, timestops)
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
    RampedVariable{true}(setpts, ramprates, holds, timestops)
end

function Base.hash(rv::RampedVariable, h::UInt)
    hash(rv.setpts, hash(rv.ramprates, hash(rv.holds, hash(rv.timestops, hash(:RampedVariable, h)))))
end

function Base.show(io::IO, rv::RampedVariable{false}) 
    return print(io, "RampedVariable($(rv.setpts))")
end
    
function Base.show(io::IO, rv::RampedVariable{true})
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

@concrete terse struct ConstPhysProp
    val
end
(cpp::ConstPhysProp)(args...) = cpp.val
Base.show(io::IO, cpp::ConstPhysProp) = print(io, "ConstPhysProp($(cpp.val))")

# function Base.show(io::IO, pp::PhysProp) 
#     return print(io, "PhysProp($(rv.setpts[1]))")
# end

"""
PrimaryDryFit: a type for storing experimental data and indicating how it should be fit.

Provided constructors:

    PrimaryDryFit(t, Tfs, Tvw, t_end)
    PrimaryDryFit(t, Tfs) = PrimaryDryFit(t, Tfs, missing, missing)
    PrimaryDryFit(t, Tfs, Tvw) = PrimaryDryFit(t, Tfs, Tvw, missing)
    PrimaryDryFit(t, Tfs, t_end::Union{Unitful.Time}, Tuple{Tt, Tt}}) where Tt  = PrimaryDryFit(t, Tfs, missing, t_end)

The use of this struct is determined in large part by the implementation of 
[`LyoPronto.obj_expT`](@ref). If a given field is not available, set it
to `missing` and things should basically work. At least `t` and `Tfs` are 
expected to always be provided.

In the end, Tfs and Tvws are each stored as a tuple of vectors, but the constructor tries to be flexible about 
allowing a single vector to be passed in place of a tuple of vectors.

`Tf_iend` and `Tvw_iend` default to `[length(Tf) for Tf in Tfs]` and `[length(Tvw) for Tvw in Tvws]`, 
respectively, with one value for each temperature series; 
they are used to dictate if a given temperature series should
be truncated sooner than the full length in the fitting procedure.
This implies that all the temperature series are valid initially at the same
time points, then stop having measured values after a different number of measurements.
If a single value is given for `Tvw`, then it is taken to be an endpoint, and `Tvw_iend` will be `missing`.

`t_end` indicates an end of drying, particularly if taken from other measurements
(e.g. from Pirani-CM convergence). If set to `missing`, it is ignored in the
objective function. If set to a tuple of two times, then in the objective function any time 
in that window is not penalized; outside that window, squared error takes over, as for the 
single time case.

Principal Cases:
- Conventional: provide only `t, Tfs`
- Conventional with Pirani ending: provide `t, Tfs, t_end`
- RF with measured vial wall: provide `t, Tfs, Tvws`, 
- RF, matching model Tvw to experimental Tf[end] without measured vial wall: provide `t, Tfs, Tvw`
"""
struct PrimaryDryFit{Tt, TT, Ti, Ttv<:AbstractVector{Tt}, TTv<:AbstractVector{TT}, 
        TTf<:Tuple{TTv, Vararg{TTv}},
        TTvw<:Union{Missing, TT, Tuple{TTv, Vararg{TTv}}},
        TTvwi<:Union{Missing, Vector{Ti}},
        Tte<:Union{Missing, Tt, Tuple{Tt, Tt}}}
    t::Ttv
    Tfs::TTf
    Tf_iend::Vector{Ti}# = [length(Tf) for Tf in Tfs]
    Tvws::TTvw# = missing
    Tvw_iend::TTvwi# = (ismissing(Tvws) ? missing : [length(Tvw) for Tvw in Tvws])
    t_end::Tte# = missing
end

# Primary constructor
function PrimaryDryFit(t, Tfs, Tvws, t_end) 
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
    PrimaryDryFit(t, Tfs, [length(Tf) for Tf in Tfs], Tvws, 
    ((ismissing(Tvws) || Tvws isa Number) ? missing : [length(Tvw) for Tvw in Tvws]),
    t_end)
end
# Convenience constructors
PrimaryDryFit(t, Tfs) = PrimaryDryFit(t, Tfs, missing, missing)
PrimaryDryFit(t, Tfs, Tvws) = PrimaryDryFit(t, Tfs, Tvws, missing)
PrimaryDryFit(t, Tfs, t_end::Union{Unitful.Time, Tuple{Tt, Tt}}) where Tt  = PrimaryDryFit(t, Tfs, missing, t_end)

function Base.:(==)(p1::PrimaryDryFit, p2::PrimaryDryFit)
    cond1 = p1.t == p2.t
    cond2 = p1.Tfs == p2.Tfs
    cond3 = p1.Tf_iend == p2.Tf_iend
    cond4 = ismissing(p1.Tvws) ? ismissing(p2.Tvws) : (p1.Tvws == p2.Tvws)
    cond5 = ismissing(p1.Tvw_iend) ? ismissing(p2.Tvw_iend) : (p1.Tvw_iend == p2.Tvw_iend)
    cond6 = ismissing(p1.t_end) ? ismissing(p2.t_end) : (p1.t_end == p2.t_end)
    return all([cond1, cond2, cond3, cond4, cond5, cond6])
end

# Add a little bit of sugar to our transforms
struct ConstWrapTV <: TransformVariables.ScalarTransform end
TransformVariables.transform(::ConstWrapTV, x) = ConstPhysProp(x)
TransformVariables.inverse(::ConstWrapTV, x) = x.value

TransformVariables.transform(t::TVScale{ConstPhysProp}, x) = ConstPhysProp(t.scale.value*x)
TransformVariables.inverse(t::TVScale{ConstPhysProp}, x) = x.value/t.scale.value