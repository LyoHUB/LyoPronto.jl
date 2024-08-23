export RpFormFit, RampedVariable

struct RpFormFit
    R0
    A1
    A2
end

function (ff::RpFormFit)(x)
    return ff.R0 + ff.A1*x/(1+ustrip(NoUnits, ff.A2*x))
end

struct RampedVariable
    setpts::AbstractVector
    ramprates::AbstractVector
    holds::AbstractVector
    timestops::AbstractVector
end

function RampedVariable(initial)
    RampedVariable([initial], [], [], [0])
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
    RampedVariable(setpts, ramprates, holds, timestops)
end

function (rv::RampedVariable)(t)
    if length(rv.timestops) == 1
        return rv.setpts[1]
    end
    im = findlast(rv.timestops .<= t)
    if im == length(rv.timestops)
        return rv.setpts[end]
    elseif iseven(im)
        return rv.setpts[im÷2+1]
    else
        ip = im+1
        return (rv.setpts[ip÷2+1] - rv.setpts[ip÷2])/(rv.timestops[ip] - rv.timestops[im])*(t - rv.timestops[im]) + rv.setpts[ip÷2]
    end
end

    
function Base.show(io::IO, rv::RampedVariable) 
    if length(rv.timestops) == 1
        return print(io, "RampedVariable($(rv.setpts[1]))")
    else 
        return print(io, "RampedVariable($(rv.setpts), $(rv.ramprates), $(rv.holds)")
    end
end