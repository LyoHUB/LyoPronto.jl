export RpFormFit, RampedVariable, ConstPhysProp

"""
A convenience type for dealing with the common functional form given to Rp and Kv.

An object `Rp = RpFormFit(A, B, C)` can be called as `Rp(x)`, which simply computes `A + B*x/(1 + C*x)`.
Likewise, `Kv = RpFormFit(Kc, Kp, Kd)` can be called as `Kv(p)` to get `Kc + Kp*p/(1 + Kd*p)`.

Be careful to pass dimensionally consistent values.
"""
struct RpFormFit{T}
    R0::T
    A1
    A2
    RpFormFit(R0, A1, A2) = RpFormFit{typeof(R0)}(R0, A1, A2)
end

function (ff::RpFormFit{T})(x) where T
    return (ff.R0 + ff.A1*x/(1+ustrip(NoUnits, ff.A2*x)))::T
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
"""
struct RampedVariable{T, vary}
    setpts::Union{T, AbstractVector{T}}
    ramprates::AbstractVector
    holds::AbstractVector
    timestops::AbstractVector
end

function RampedVariable(initial)
    RampedVariable{typeof(initial), false}([initial], [], [], [0])
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
    RampedVariable{eltype(setpts), true}(setpts, [ramprate], [], timestops)
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
    for (i, ramp) in enumerate(rest)
        timestops[2i+1] = timestops[2i] + holds[i]
        timestops[2i+2] = timestops[2i+1] + (setpts[i+2]-setpts[i+1])/ramp
    end
    RampedVariable{eltype(setpts), true}(setpts, ramprates, holds, timestops)
end

function (rv::RampedVariable{T, false})(t) where T
    return rv.setpts[1]::T
end
function (rv::RampedVariable{T, true})(t) where T
    im = findlast(rv.timestops .<= t)
    if im == length(rv.timestops)
        return rv.setpts[end]
    elseif iseven(im)
        return rv.setpts[im÷2+1]
    else
        ip = im+1
        return ((rv.setpts[ip÷2+1] - rv.setpts[ip÷2])/(rv.timestops[ip] - rv.timestops[im])*(t - rv.timestops[im]) + rv.setpts[ip÷2])::T
    end
end

    
function Base.show(io::IO, rv::RampedVariable) 
    if length(rv.timestops) == 1
        return print(io, "RampedVariable($(rv.setpts[1]))")
    else 
        return print(io, "RampedVariable($(rv.setpts), $(rv.ramprates), $(rv.holds)")
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