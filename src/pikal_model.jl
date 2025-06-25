export lyo_1d!, lyo_1d_dae_f, subflux_Tsub, calc_psub
export end_drying_callback
export ParamObjPikal
export calc_u0, get_tstops
export RpEstimator, calc_hRp_T

@doc raw"""
    end_cond(u, t, integ)

Compute the end condition for primary drying (that `mf` or `hf` approaches zero).
"""
end_cond(u, t, integ) = u[1] - 1e-10 # When reaches 1e-10, is basically zero
@doc raw"""
A callback for use in simulating either the Pikal or RF model.

Terminates the time integration when [`end_cond`](@ref) evaluates to `true`.
"""
const end_drying_callback = ContinuousCallback(end_cond, terminate!, save_positions=(true, false))

# -------------------------------------------
# Incorporate the nonlinear algebraic part in a DAE formulation.
# This has the advantage that, afterward, temperatures can be cheaply interpolated by builtin solutions

@doc raw"""
    lyo_1d_dae!(du, u, params, t)

Internal implementation of the Pikal model.
See [`lyo_1d_dae_f`](@ref) for the wrapped version, which is more fully documented.
"""
function lyo_1d_dae!(du, u, params, t)

    if params isa ParamObj
        (; Rp, hf0, csolid, ρsolution,
        Kshf, Av, Ap, pch, Tsh) = params
    else
        Rp, hf0, csolid, ρsolution = params[1]
        Kshf, Av, Ap = params[2]
        pch, Tsh = params[3]
    end
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    if Tf < 0.0u"K"
        du .= NaN
        return nothing
    end

    hd = hf0 - hf
    pchl = pch(td)
    Qshf = Av*Kshf(pchl)*(Tsh(td) - Tf)
    Tsub = Tf - Qshf/k_ice/Ap*hf
    delta_p = calc_psub(Tsub)-pch(td)
    dmdt = - Ap*(delta_p)/Rp(hd)
    Qsub = uconvert(u"W", dmdt*ΔHsub)

    dhf_dt = min(0.0u"cm/hr", dmdt/(ρsolution-csolid)/Ap |> u"cm/hr") # Cap dhf_dt at 0: no desublimation

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Qsub + Qshf)
    return nothing
end

const lyo_1d_mm = Diagonal([1.0, 0.0])

@doc raw"""
    lyo_1d_dae_f = ODEFunction(lyo_1d_dae!, mass_matrix=lyo_1d_mm)

Compute the right hand side function for the Pikal model.

The DAE system which is the Pikal model (1 ODE, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_1d_dae!`](@ref)

The initial conditions `u0 = [h_f, Tf]` should be unitless, but are internally assigned to be in `[cm, K]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr]`.

`params` is a `ParamObjPikal`, which can be constructed with the following form (helping with readability):
```
params = ParamObjPikal((
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh) ,
))
```
where those listed following are callables returning `Quantity`s, and the rest are `Quantity`s.
See [`RpFormFit`](@ref LyoPronto.RpFormFit) and [`RampedVariable`](@ref LyoPronto.RampedVariable) for convenience types that can help with the callables.
- `Rp(x)` with `x` a length returns mass transfer resistance (as a Unitful quantity)
- `Kshf(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)` return shelf temperature and chamber pressure respectively at time `t`.
"""
const lyo_1d_dae_f = ODEFunction{true, SciMLBase.AutoSpecialize}(lyo_1d_dae!, mass_matrix=lyo_1d_mm)


# ```
# params = (
#     (Rp, hf0, csolid, ρsolution),
#     (Kshf, Av, Ap),
#     (pch, Tsh) ,
# )
# ```
@concrete terse struct ParamObjPikal <: ParamObj
    Rp
    hf0
    csolid
    ρsolution
    Kshf
    Av
    Ap
    pch
    Tsh
end

@doc """
    $(TYPEDEF)

The `ParamObjPikal` type is a container for the parameters used in the Pikal model.
"""
ParamObjPikal

# This constructor takes the legacy tuple of tuples form I used and unpacks it
function ParamObjPikal(tuptup) 
    return ParamObjPikal(tuptup[1]..., tuptup[2]..., tuptup[3]...)
end

# This constructor catches if Rp is passed as a tuple or NamedTuple and constructs an appropriate object
function ParamObjPikal(Rp::Union{NamedTuple, Tuple},
    hf0, csolid, ρsolution, Kshf, Av, Ap, pch, Tsh)
    return ParamObjPikal(RpFormFit(Rp...), hf0, csolid, ρsolution, Kshf, Av, Ap, pch, Tsh)
end

function Base.getindex(p::ParamObjPikal, i::Int)
    if i == 1
        return (p.Rp, p.hf0, p.csolid, p.ρsolution)
    elseif i == 2
        return (p.Kshf, p.Av, p.Ap)
    elseif i == 3
        return (p.pch, p.Tsh)
    else
        error(BoundsError, "Attempt to access LyoPronto.ParamsObjPikal at index $i. Only indices 1 to 3 allowed")
    end
end
Base.size(::ParamObjPikal) = (3,)
Base.length(::ParamObjPikal) = 3

function calc_u0(po::ParamObjPikal)
    # return ustrip(u"cm", u"K"),[po.hf0, po.Tsh(0u"s")])
    return [ustrip(u"cm", po.hf0), ustrip(u"K", float(po.Tsh(0u"s")))]
end
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
function get_tstops(po::ParamObjPikal)
    get_tstops((po.Tsh, po.pch))
end

function get_t0(Tsh, pch)
    if calc_psub(Tsh(0u"s")) < pch(0u"s")
        t0 = find_zero(t -> ustrip(u"Pa", calc_psub(Tsh(t*u"hr")) - pch(t*u"hr")), (0.0, ustrip(u"hr", Tsh.timestops[end])))
        return t0 * 1.001 # Go jslightly after the zero, to ensure stability
    else
        return 0.0
    end
end
get_t0(po::ParamObjPikal) = get_t0(po.Tsh, po.pch)

function ODEProblem(po::ParamObjPikal; u0=calc_u0(po), tspan=(0.0, 1000.0))
    tstops = get_tstops(po)
    t0 = get_t0(po)
    @reset tspan[1] = t0 
    return ODEProblem(lyo_1d_dae_f, u0, tspan, po; 
        tstops = tstops, callback=end_drying_callback)
end

# -------------------------------------------
# Alternative approach, a la original LyoPRONTO
# Embed the nonlinear solve inside the function call for a plain ODE formulation.

function compute_T_pseudosteady(Pchl, Rpl, Kvl, Tshl, Ap, Av, hd)
    function nl_func(Tp_nd)
        Tp = Tp_nd*u"K"
        Qs = Av*Kvl*(Tshl - Tp)
        Tsub = Tp - Qs/Ap/k_ice*hd
        dmdt = - max(Ap*(calc_psub(Tsub)-Pchl)/Rpl, 0u"kg/s")
        Qsub = dmdt*ΔHsub
        return ustrip(u"W", Qsub + Qs)
    end
    Tp = find_zero(nl_func, 250) *u"K"
end

function lyo_1d!(du, u, params, t)
    Rp, hf0, csolid, ρsolution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    hd = hf0 - hf

    Tp = compute_T_pseudosteady(pch(td), Rp(hd), Kshf(pch(td)), Tsh(td), Ap, Av, hd)
    dmdt = - Ap*(calc_psub(Tp)-pch(td))/Rp(hd)

    dhf_dt = min(0u"cm/s", dmdt/Ap/(ρsolution - csolid))

    # Q = uconvert(u"W/m^2", Kshf(pch(td))*(Tsh(td)-Tp))
    # flux = uconvert(u"kg/hr/m^2", -dmdt/Ap)

    du[1] = min(0, ustrip(u"cm/hr", dhf_dt))
end

function subflux_Tsub(u, params, t)
    Rp, hf0, csolid, ρsolution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 

    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"

    hd = hf0 - hf
    Qshf = Av*Kshf(pch(td))*(Tsh(td) - Tf)
    Tsub = Tf - Qshf/k_ice/Ap*hf
    dmdt = - max(Ap*(calc_psub(Tsub)-pch(td))/Rp(hd), 0u"kg/s")
    # Qsub = dmdt*ΔHsub
    return dmdt, Tsub
end

function subflux_Tsub(sol::T, t) where T<:ODESolution
    subflux_Tsub(sol(t), sol.p, t)
end

# -----------------
# Directly estimate Rp from time series

struct RpEstimator{plural}
    po::ParamObjPikal
    pdf::PrimaryDryFit
    Tf_interp
end

function RpEstimator(po::ParamObjPikal, pdf::PrimaryDryFit)
    if length(pdf.Tf_iend) == 1
        return RpEstimator{false}(po, pdf, LinearInterpolation(pdf.Tfs[1], pdf.t[begin:pdf.Tf_iend[1]]))
    end
    Tf_interp = [LinearInterpolation(pdf.Tfs[i], pdf.t[begin:i_end], extrapolation=ExtrapolationType.Constant) for (i, i_end) in enumerate(pdf.Tf_iend)]
    return RpEstimator{true}(po, pdf, Tf_interp)
end

function Base.show(io::IO, re::RpEstimator{plural}) where plural
    return print(io, "RpEstimator{$plural}(...)")
end

function Base.getindex(re::RpEstimator{true}, i)
    return RpEstimator{false}(re.po, re.pdf, re.Tf_interp[i])
end
Base.length(re::RpEstimator{false}) = length(re.Tf_interp.t)


function dae_Rp!(du, u, p, tn)
    t = tn*u"hr"
    hd = u[1]*u"cm"
    Rpg = u[2]*u"cm^2*Torr*hr/g"

    (;po, Tf_interp) = p
    (;hf0, csolid, ρsolution, 
    Kshf, Av, Ap, 
    pch, Tsh) = po

    Tf = Tf_interp(t)


    Q = Kshf(pch(t))*Av*(Tsh(t) - Tf)
    Tsub = Tf - Q/Ap/LyoPronto.k_ice * (hf0-hd)
    md = Q/LyoPronto.ΔH
    Rp = Ap*(calc_psub(Tsub)-pch(t))/md |> u"cm^2*Torr*hr/g"
    if Q <= 0.0u"W" || Rp <= 0.0u"m/s" || isnan(Rp)
        du[1] = du[2] = 0.0
        return
    end

    du[1] = ustrip(u"cm/hr", md/(ρsolution-csolid)/Ap)
    du[2] = u[2] - ustrip(u"cm^2*Torr*hr/g", Rp)
    return
end
const dae_Rpf = ODEFunction(dae_Rp!, mass_matrix=Diagonal([1.0, 0]))

function get_t0(re::RpEstimator{false})
    (; Tf_interp, po) = re
    t0 = 0.0
    if calc_psub(Tf_interp(0u"hr")) < po.pch(0u"hr")
        t0_psub = find_zero((0.0, ustrip(u"hr", Tf_interp.t[end]))) do t
            Tf = Tf_interp(t*u"hr")
            Tsh = po.Tsh(t*u"hr")
            pch = po.pch(t*u"hr")
            Q = po.Kshf(pch)*po.Av*(Tsh - Tf)
            Tsub = Tf - Q/po.Ap/LyoPronto.k_ice*po.hf0
            ustrip(u"Pa", calc_psub(Tsub)-pch)
        end
        t0 = max(t0, t0_psub*1.01)
    end
    if Tf_interp(0u"hr") > po.Tsh(0u"hr")
        t0_q = find_zero((0.0, ustrip(u"hr", Tf_interp.t[end]))) do t
            ustrip(u"K", Tf_interp(t*u"hr") - po.Tsh(t*u"hr"))
        end
        t0 = max(t0, t0_q*1.01)
    end
    return t0
end

ODEProblem(::RpEstimator{true}) = error("Cannot create ODEProblem for multiple Tf at once. Index into the RpEstimator to choose a Tf series.")
function ODEProblem(re::RpEstimator{false}; u0=[0.0,0], tspan=(get_t0(re), ustrip(u"hr", re.Tf_interp.t[end])))
    return ODEProblem(dae_Rpf, u0, tspan, re; tstops=ustrip.(u"hr", re.Tf_interp.t))
end

function calc_hRp_T(po::ParamObjPikal, pdf::PrimaryDryFit; i=nothing)
    re = RpEstimator(po, pdf)
    if re isa RpEstimator{true}
        if !isnothing(i)
            prob = ODEProblem(re[i])
        else
            @warn "Index needed for multiple Tf. Taken as 1 by default" i 
            prob = ODEProblem(re[1])
        end
    else
        !isnothing(i) && @warn "Index passed but not needed" i 
        prob = ODEProblem(re)
    end
    sol = solve(prob, Rosenbrock23(), saveat=ustrip.(u"hr", pdf.t))
    hd, Rp = sol[1,:]*u"cm", sol[2,:]*u"cm^2*Torr*hr/g"
    # If there are multiple zeros at the start, trim them off
    hdi = findlast(hd .== 0.0u"cm")
    Rpi = findlast(Rp .== 0.0u"m/s")
    if !isnothing(hdi) && !isnothing(Rpi) 
        istart = max(min(hdi-1, Rpi-1), 1)
        hd = hd[istart:end]
        Rp = Rp[istart:end]
    end
    return hd, Rp
end
