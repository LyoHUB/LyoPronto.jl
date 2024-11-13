export lyo_1d!, lyo_1d_dae_f, subflux_Tsub, calc_psub
export end_drying_callback

@doc raw"""
    end_cond(u, t, integ)

Compute the end condition for primary drying (that `mf` or `hf` approaches zero).
"""
end_cond(u, t, integ) = u[1] - 1e-10 # When reaches 1e-10, is basically zero
@doc raw"""
A callback for use in simulating either the Pikal or RF model.

Terminates the time integration when [`end_cond`](@ref) evaluates to `true`.
"""
const end_drying_callback = ContinuousCallback(end_cond, terminate!)

# -------------------------------------------
# Incorporate the nonlinear algebraic part in a DAE formulation.
# This has the advantage that, afterward, temperatures can be cheaply interpolated by builtin solutions

@doc raw"""
    lyo_1d_dae!(du, u, params, t)

Internal implementation of the Pikal model.
See [`lyo_1d_dae_f`](@ref) for the wrapped version, which is more fully documented.
"""
function lyo_1d_dae!(du, u, params, t)

    Rp, hf0, c_solid, ρ_solution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 
    # Deliberately ignore other parameters that may be used for RF model
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    if Tf < 0u"K"
        return [NaN, NaN]
    end

    hd = hf0 - hf
    pchl = pch(td)
    Qshf = Av*Kshf(pchl)*(Tsh(td) - Tf)
    Tsub = Tf - Qshf/k_ice/Ap*hf
    dmdt = - Ap*(calc_psub(Tsub)-pch(td))/Rp(hd)
    Qsub = uconvert(u"W", dmdt*ΔHsub)

    dhf_dt = min(0u"cm/hr", dmdt/(ρ_solution-c_solid)/Ap) # Cap dhf_dt at 0: no desublimation

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Qsub + Qshf)
end

const lyo_1d_mm = [1.0 0.0; 0.0 0.0]

@doc raw"""
    lyo_1d_dae_f = ODEFunction(lyo_1d_dae!, mass_matrix=lyo_1d_mm)

Compute the right hand side function for the Pikal model.

The DAE system which is the Pikal model (1 ODE, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_1d_dae!`](@ref)

The initial conditions `u0 = [h_f, Tf]` should be unitless, but are internally assigned to be in `[cm, K]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr]`.

`params` has the form:
```
params = (
    (Rp, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh) ,
)
```
where those listed following are callables returning `Quantity`s, and the rest are `Quantity`s.
See [`RpFormFit`](@ref LyoPronto.RpFormFit) and [`RampedVariable`](@ref LyoPronto.RampedVariable) for convenience types that can help with the callables.
- `Rp(x)` with `x` a length returns mass transfer resistance (as a Unitful quantity)
- `Kshf(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)` return shelf temperature and chamber pressure respectively at time `t`.
"""
lyo_1d_dae_f = ODEFunction(lyo_1d_dae!, mass_matrix=lyo_1d_mm)

# -------------------------------------------
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
    Rp, hf0, c_solid, ρ_solution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    hd = hf0 - hf

    Tp = compute_T_pseudosteady(pch(td), Rp(hd), Kshf(pch(td)), Tsh(td), Ap, Av, hd)
    dmdt = - Ap*(calc_psub(Tp)-pch(td))/Rp(hd)

    dhf_dt = min(0u"cm/s", dmdt/Ap/(ρ_solution - c_solid))

    # Q = uconvert(u"W/m^2", Kshf(pch(td))*(Tsh(td)-Tp))
    # flux = uconvert(u"kg/hr/m^2", -dmdt/Ap)

    du[1] = min(0, ustrip(u"cm/hr", dhf_dt))
end

function subflux_Tsub(u, params, t)
    Rp, hf0, c_solid, ρ_solution = params[1]
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