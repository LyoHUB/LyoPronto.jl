export lumped_cap_rf!
export ParamObjRF

function shapefac(Bi)
    charfunc(x) = -x*besselj1(x) + Bi*besselj0(x)
    λm = zeros(200)
    for i in eachindex(λm)
        if i == 1
            λm[i] = find_zero(charfunc, π/2)
        else
            λm[i] = find_zero(charfunc, λm[i-1]+π)
        end
    end
    Cm = @. 2/λm*besselj1(λm)/(besselj0(λm)^2 + besselj1(λm)^2)
    integ = sum(Cm .* besselj1.(λm) .* tanh.(λm))
end

const Bi_samp = 10.0 .^range(-2, 5, length=71)
const S_samp = shapefac.(Bi_samp)
const S_interp = LinearInterpolation(S_samp, Bi_samp, extrapolation=ExtrapolationType.Linear)

"""
    $(SIGNATURES)

Compute the right-hand-side function for the ODEs making up the lumped-capacitance microwave-assisted model.

The optional argument `qret` defaults to `Val(false)`; if set to `Val(true)`, the function returns
`[Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw]` with `Q_...` as Unitful quantities in watts. 
The extra results are helpful in investigating the significance of the various heat transfer 
modes, but are not necessary in the ODE integration.

`du` refers to `[dmf/dt, dTf/dt, dTvw/dt]`, with `u = [mf, Tf, Tvw]`.
`u` is taken without units but assumed to have the units of `[g, K, K]` (which is internally added).
`tn` is assumed to be in hours (internally added), so `dudt` is returned with assumed units `[g/hr, K/hr, K/hr]` to be consistent.

It is recommended to use the `ParamObjRF` type to hold the parameters, since it allows some
more convenient access to the parameters, but they can be given in the form of a 
tuple-of-tuples:
```
params = (   
    (Rp, hf0, cSolid, ρsolution),
    (Kshf_f, Av, Ap),
    (pch, Tsh, P_per_vial),
    (mf0, cpf, mv, cpv, Arad),
    (f_RF, eppf, eppvw),
    (Kvwf, Bf, Bvw, alpha),
)
```
This tuple-of-tuples structure is also used in an extra constructor for the `ParamObjRF` type.

These should all be Unitful quantities with appropriate dimensions, with some exceptions which are callables returning quantities.
See [`RpFormFit`](@ref) and [`RampedVariable`](@ref) for convenience types that can help with these cases.
- `Rp(x)` with `x` a length returns mass transfer resistance (as a Unitful quantity)
- `Kshf_f(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)`, `P_per_vial(t)` return shelf temperature, chamber pressure, and microwave power respectively at time `t`.

- `Arad` and `alpha` are used only in prior versions of the model, and can be left out.
This is my updated version of the model.
LC3: Q_shw evaluated with Kshf; shape factor included; α=0
"""
function lumped_cap_rf!(du, u, params, tn, qret = Val(false))
    # Unpack all the parameters
    if params isa ParamObjRF
        (;Rp, hf0, csolid, ρsolution,
        Kshf, Av, Ap,
        pch, Tsh, P_per_vial, 
        mf0, cpf, mv, cpv,
        f_RF, eppf, eppvw,
        Kvwf, Bf, Bvw) = params
    else
        Rp, hf0, csolid, ρsolution = params[1]
        Kshf, Av, Ap, = params[2]
        pch, Tsh, P_per_vial = params[3] 
        mf0, cpf, mv, cpv = params[4]
        f_RF, eppf, eppvw = params[5]
        Kvwf, Bf, Bvw = params[6]
    end
    # Dimensionalize the state variables
    t = tn*u"hr" 
    m_f = u[1]*u"g"
    T_f = u[2]*u"K"
    T_vw = u[3]*u"K"
    # Compute some properties
    porosity = (ρsolution - csolid)/ρsolution
    k_dry = k_sucrose*(1-porosity)
    V_vial = mv / rho_glass
    # Do some geometry
    rad = sqrt(Ap/π)
    h_f = m_f/mf0 * hf0 
    h_d = hf0 - h_f
    # Heat transfer from shelf
    Kshft = Kshf(pch(t))
    Q_shf = Kshft*Ap*(Tsh(t)-T_f) 
    Q_shw = Kshft*(Av-Ap)*(Tsh(t)-T_vw)
    # Evaluate mass flow; positive means drying is progressing. Not forced to be positive
    mflow = Ap/Rp(h_d)*(calc_psub(T_f) - pch(t)) # g/s
    Q_sub = mflow*ΔHsub # Sublimation
    # Evaluate heat transfer from wall
    Bi = uconvert(NoUnits, Kvwf*rad/k_dry)
    Q_vwf = 2π*(Kvwf*rad*h_f + k_dry*(hf0-h_f)*S_interp(Bi)) * (T_vw-T_f)
    # Volumetric heating
    Qppp_RF_f  = 2*pi*f_RF*e_0*eppf(T_f, f_RF)*P_per_vial(t)*Bf # W / m^3
    Qppp_RF_vw = 2*pi*f_RF*e_0*eppvw*P_per_vial(t)*Bvw # W / m^3
    Q_RF_f = Qppp_RF_f*Ap*h_f # W
    Q_RF_vw = Qppp_RF_vw*V_vial # W
    # Check that total volumetric heating is less than input power
    if Q_RF_f + Q_RF_vw > P_per_vial(t) && t == 0u"hr"
        @warn "Energy balance of EM terms not satisfied." uconvert(u"W", Q_RF_f) uconvert(u"W", Q_RF_vw) P_per_vial(t)
    end
    # Evaluate derivatives
    # Desublimation is not allowed here: if we clamp mflow itself, then the DAE is unstable
    dm_f = min(0.0u"kg/s", -mflow/porosity)
    dT_f =  (Q_shf+Q_vwf+Q_RF_f -Q_sub) / (m_f*cpf) - T_f*dm_f/m_f
    dT_vw = (Q_shw-Q_vwf+Q_RF_vw) / (mv*cpv)

    du[1] = ustrip(u"g/hr", dm_f)
    du[2] = ustrip(u"K/hr", dT_f)
    du[3] = ustrip(u"K/hr", dT_vw)
    # Strip units from derivatives; return all heat transfer terms
    if qret isa Val{true}
        return uconvert.(u"W", [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw])
    else
        return nothing
    end
end

# ```
# params = (   
#     (Rp, hf0, csolid, ρsolution),
#     (Kshf, Av, Ap),
#     (pch, Tsh, P_per_vial),
#     (mf0, cpf, mv, cpv, Arad),
#     (f_RF, eppf, eppvw),
#     (Kvwf, Bf, Bvw, alpha),
# )
# ```

@concrete terse struct ParamObjRF <: ParamObj
    Rp
    hf0
    csolid
    ρsolution
    Kshf
    Av
    Ap
    pch
    Tsh
    P_per_vial
    mf0
    cpf
    mv
    cpv
    Arad
    f_RF
    eppf
    eppvw
    Kvwf
    Bf
    Bvw
    alpha
end

@doc """
    $(TYPEDEF)

The `ParamObjRF` type is a container for the parameters used in the RF model.
"""
ParamObjRF

function ParamObjRF(tuptup) 
    if (length(tuptup[4]) == 4 && length(tuptup[6]) == 3)
        return ParamObjRF(tuptup[1]..., tuptup[2]...,
                    tuptup[3]..., tuptup[4]..., missing,
                    tuptup[5]..., tuptup[6]..., missing,)
    elseif length(tuptup[6]) == 3
        return ParamObjRF(tuptup[1]..., tuptup[2]...,
                    tuptup[3]..., tuptup[4]...,
                    tuptup[5]..., tuptup[6]..., missing,)
    else
        return ParamObjRF(tuptup[1]..., tuptup[2]...,
                    tuptup[3]..., tuptup[4]...,
                    tuptup[5]..., tuptup[6]...,)
    end
end
Base.size(po::ParamObjRF) = (6,)

function Base.getindex(po::ParamObjRF, i)
    if i == 1
        return (po.Rp, po.hf0, po.csolid, po.ρsolution)
    elseif i == 2
        return (po.Kshf, po.Av, po.Ap)
    elseif i==3 
        return (po.pch, po.Tsh, po.P_per_vial)
    elseif i == 4
        return (po.mf0, po.cpf, po.mv, po.cpv, po.Arad)
    elseif i == 5
        return (po.f_RF, po.eppf, po.eppvw)
    elseif i == 6
        return (po.Kvwf, po.Bf, po.Bvw, po.alpha)
    else
        error(BoundsError, "Attempt to access LyoPronto.ParamsObjRF at index $i. Only indices 1 to 6 allowed")
    end
end

function calc_u0(po::ParamObjRF)
    Tsh0_nd = ustrip(u"K", float(po.Tsh(0u"s")))
    # return ustrip.((u"g", u"K", u"K"), (po.mf0, Tsh0, Tsh0))
    return [ustrip(u"g", po.mf0), Tsh0_nd, Tsh0_nd]
end
function get_tstops(po::ParamObjRF)
    get_tstops((po.Tsh, po.pch, po.P_per_vial))
end

function ODEProblem(po::ParamObjRF; u0 = calc_u0(po), tspan=(0.0, 400.0))
    tstops = get_tstops(po)
    return ODEProblem{true, SciMLBase.FullSpecialize}(lumped_cap_rf!, u0, tspan, po; 
        tstops = tstops, callback=end_drying_callback)
end