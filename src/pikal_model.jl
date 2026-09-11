
# -------------------------------------------
# Incorporate the nonlinear algebraic part in a DAE formulation.
# This has the advantage that, afterward, temperatures can be cheaply interpolated by builtin solutions

const PIKAL_PARAM_DOC = """
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


"""
    $(SIGNATURES)

With the Pikal model, compute model quantities at `u=[hf, Tf]`, time `t`, and conditions `po`.

`po` should be a [`ParamObjPikal`](@ref); the return will be a named tuple with Unitful quantities in the following fields:
- `md`: mass flow rate (g/hr)
- `Q_shf`: heat transfer from shelf to product (W)

This allows assessment of the model's outputs without needing to rewrite the model equations.
"""
@inline function calc_md_Q(u, po, t)

    (; Rp, hf0, csolid, ρsolution,
    Kshf, Av, Ap, pch, Tsh) = po

    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    hd = hf0 - hf
    # escape hatch: unphysical temperatures or Rp too small can cause DAE convergence failure
    if Tf < 0.0u"K" || Rp(hd) < 1e-4u"hr*cm^2*Torr/g"
        return (; md=NaN*u"kg/s", Q_shf=NaN*u"W")
    end

    pchl = pch(td)
    Q_shf = Av*Kshf(pchl)*(Tsh(td) - Tf) |> u"W"
    Tsub = Tf - Q_shf/k_ice/Ap*hf
    delta_p = calc_psub(Tsub)-pch(td)
    md = - Ap*(delta_p)/Rp(hd) |> u"g/hr"
    return (; md, Q_shf)
end

@doc raw"""
    lyo_1d_dae!(du, u, params, t)

Internal implementation of the Pikal model.
See [`lyo_1d_dae_f`](@ref) for the wrapped version, which is more fully documented.
"""
function lyo_1d_dae!(du, u, params, t)
    
    # Need a handful of parameters in this function.
    if params isa ParamObj
        (; csolid, ρsolution, Ap) = params
    else
        csolid, ρsolution = params[1][3:4]
        Ap = params[2][3]
    end
    # This logic is carried out in a separate function,
    # so that it can be reused after the fact for computing mass flow.
    (;md, Q_shf) = calc_md_Q(u, params, t)
    dmdt = md
    if isnan(dmdt)
        du .= NaN
        return nothing
    end
    Q_sub = uconvert(u"W", dmdt*ΔHsub)

    dhf_dt = min(0.0u"cm/hr", dmdt/(ρsolution-csolid)/Ap |> u"cm/hr") # Cap dhf_dt at 0: no desublimation

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Q_sub + Q_shf)
    return nothing
end

const lyo_1d_mm = Diagonal([1.0, 0.0])

"""
    lyo_1d_dae_f = ODEFunction(lyo_1d_dae!, mass_matrix=Diagonal([1.0, 0.0]))

Compute the right hand side function for the Pikal model.

The DAE system which is the Pikal model (1 ODE, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_1d_dae!`](@ref) and [`calc_md_Q`](@ref).

The initial conditions `u0 = [h_f, Tf]` should be unitless, but are internally assigned to be in `[cm, K]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr]`.

$(PIKAL_PARAM_DOC)
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

$(PIKAL_PARAM_DOC)
"""
ParamObjPikal

# This constructor takes the legacy tuple of tuples form I used and unpacks it
function ParamObjPikal(tuptup) 
    return ParamObjPikal(tuptup[1]..., tuptup[2]..., tuptup[3]...)
end

# -------------------------------------------
# Define how a ParamObjPikal maps to an ODEProblem

# Define how u0, tstops, and initial time should be calculated for a given ParamObjPikal
function calc_u0(po::ParamObjPikal)
    return [ustrip(u"cm", po.hf0), ustrip(u"K", float(po.Tsh(0u"s")))]
end
function get_tstops(po::ParamObjPikal)
    get_tstops((po.Tsh, po.pch))
end
get_t0(po::ParamObjPikal) = get_t0(po.Tsh, po.pch)

function ODEProblem(po::ParamObjPikal; u0=calc_u0(po), tspan=(0.0, 1000.0))
    tstops = get_tstops(po)
    t0 = get_t0(po)
    @reset tspan[1] = t0 
    return ODEProblem{true, SciMLBase.FullSpecialize}(lyo_1d_dae_f, u0, tspan, po; 
        tstops = tstops, callback=end_drying_callback, initializealg=BrownFullBasicInit(),
        dt=0.1)
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
    Rp = uconvert(u"cm^2*Torr*hr/g", Ap*(calc_psub(Tsub)-pch(t))/md)
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
    tr = Tf_interp.t
    t0 = tr[1]
    for i in 1:(length(tr)÷2) # check first half of data points
        t = tr[i]
        Tf = Tf_interp(t)
        Tsh = po.Tsh(t)
        pch = po.pch(t)
        Q = po.Kshf(pch)*po.Av*(Tsh - Tf)
        if Q < 0u"W" # Negative heat transfer
            t0 = tr[i+1]
            continue
        end
        Tsub = Tf - Q/po.Ap/LyoPronto.k_ice*po.hf0
        if calc_psub(Tsub) < pch # Negative mass transfer
            t0 = tr[i+1]
        end
    end
    return ustrip(u"hr", t0)
end

ODEProblem(::RpEstimator{true}) = error("Cannot create ODEProblem for multiple Tf at once. Index into the RpEstimator to choose a Tf series.")
function ODEProblem(re::RpEstimator{false}; u0=[0.0,0], tspan=(get_t0(re), ustrip(u"hr", re.Tf_interp.t[end])))
    return ODEProblem(dae_Rpf, u0, tspan, re; tstops=ustrip.(u"hr", re.Tf_interp.t), initializealg=BrownFullBasicInit())
end

"""
    $(SIGNATURES)

For experimental conditions given by a `po` and experimental data given by `pdf`, 
compute the effective \$R_p\$ and dry layer height \$h_d\$ over time.

Since `po` is a a `ParamObjPikal`, you will need to construct that object--the value of \$R_p\$ 
will not be used here, so set it to any dummy value. 

If `pdf` has multiple temperature series, pass an index `i` to select which series to use. 
Otherwise, the first series will be used.
"""
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
    sol = solve(prob, odealg_chunk2, saveat=ustrip.(u"hr", pdf.t))
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
