
# ---------------------------------
# Callback for end of drying
"""
    end_cond(u, t, integ)

Compute the end condition for primary drying (that `mf` or `hf` approaches zero).
"""
end_cond(u, t, integ) = u[1] - 1e-10 # When reaches 1e-10, is basically zero
"""
A callback for use in simulating either the Pikal or RF model.

Terminates the time integration when [`end_cond`](@ref) evaluates to `true`.
"""
const end_drying_callback = ContinuousCallback(end_cond, terminate!, save_positions=(true, false))

# -------------------------------------------
# Functional form for Rp, as a callable object

@concrete terse struct RpFormFit
    R0
    A1
    A2
end
RpFormFit(;R0, A1, A2) = RpFormFit(R0, A1, A2)
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

# ---------------------------------------------------------------
# RampedVariable:

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

# ----------------------------------------------
# Helpers for getting time stops and initial time

extract_ts(rv::RampedVariable{true, T1, T2, T3, T4}; un=u"hr") where {T1, T2, T3, T4} = ustrip.(un, float.(rv.timestops))
extract_ts(rv::RampedVariable{false, T1, T2, T3, T4}; un=u"hr") where {T1, T2, T3, T4} = [0.0]
extract_ts(interp::DataInterpolations.AbstractInterpolation; un=u"hr") = ustrip.(un, float.(interp.t))
extract_ts(a::Any) = [0.0]
function get_tstops(controls::Tuple)
    # tstops = [0.0]
    tstops = mapreduce(extract_ts, vcat, controls)
    # tstops = vcat(tstops, newstops)
    sort!(tstops); unique!(tstops)
    return tstops
end

function get_t0(Tsh, pch)
    if calc_psub(Tsh(0u"s")) < pch(0u"s")
        t0 = find_zero(t -> ustrip(u"Pa", calc_psub(Tsh(t*u"hr")) - pch(t*u"hr")), (0.0, ustrip(u"hr", Tsh.timestops[end])))
        return t0 * 1.001 # Go slightly after the zero, to ensure stability
    else
        return 0.0
    end
end

# --------------------------------------
# Callables that wrap a single constant value
@concrete terse struct ConstPhysProp
    val
end
(cpp::ConstPhysProp)(args...) = cpp.val
Base.show(io::IO, cpp::ConstPhysProp) = print(io, "ConstPhysProp($(cpp.val))")
