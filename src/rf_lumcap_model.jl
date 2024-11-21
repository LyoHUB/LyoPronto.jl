export lumped_cap_rf, lumped_cap_rf_LC1
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
const S_interp = LinearInterpolation(S_samp, Bi_samp, extrapolate=true)


@doc raw"""
    lumped_cap_rf(u, params, tn)

Compute the right-hand-side function for the ODEs making up the lumped-capacitance microwave-assisted model.

Specifically, this is `[dmf/dt, dTf/dt, dTvw/dt]` given `u = [mf, Tf, Tvw]`.
`u` is taken without units but assumed to have the units of `[g, K, K]` (which is internally added).
`tn` is assumed to be in hours (internally added), so `dudt` is returned with assumed units `[g/hr, K/hr, K/hr]` to be consistent.

The full set of necessary parameters is given in the form of a tuple-of-tuples:
```
params = (   
    (Rp, h_f0, cSolid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (m_f0, cp_f, m_v, cp_v, A_rad),
    (f_RF, epp_f, epp_vw),
    (K_vwf, B_f, B_vw, alpha),
)
```

These should all be Unitful quantities with appropriate dimensions, with some exceptions which are callables returning quantities.
See [`RpFormFit`](@ref) and [`RampedVariable`](@ref) for convenience types that can help with these cases.
- `Rp(x)` with `x` a length returns mass transfer resistance (as a Unitful quantity)
- `K_shf_f(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)`, `P_per_vial(t)` return shelf temperature, chamber pressure, and microwave power respectively at time `t`.

- `A_rad` and `alpha` are used only in the LC2 and LC3 versions of the model, and can be left out.

For implementation details, see [`lumped_cap_rf_LC1`](@ref).
"""
function lumped_cap_rf(u, params, tn)
    lumped_cap_rf_LC1(u, params, tn)[1]
end


@doc raw"""
    lumped_cap_rf_LC1(u, params, tn)

This does the work for [`lumped_cap_rf`](@ref), but returns `dudt,  [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw]` with `Q_...` as Unitful quantities in watts. 
The extra results are helpful in investigating the significance of the various heat transfer modes,
but are not necessary in the ODE integration.

LC1: Q_shw evaluated with K_shf; shape factor included; α=0
My preferred version.
"""
function lumped_cap_rf_LC1(u, params, tn)
    # Unpack all the parameters
    Rp, h_f0, cSolid, ρ_solution = params[1]
    K_shf_f, A_v, A_p, = params[2]
    pch, Tsh, P_per_vial = params[3] 
    m_f0, cp_f, m_v, cp_v = params[4]
    f_RF, epp_f, epp_vw = params[5]
    K_vwf, B_f, B_vw = params[6]
    # Dimensionalize the state variables
    t = tn*u"hr" 
    m_f = u[1]*u"g"
    T_f = u[2]*u"K"
    T_vw = u[3]*u"K"
    # Compute some properties
    porosity = (ρ_solution - cSolid)/ρ_solution
    k_dry = k_sucrose*(1-porosity)
    V_vial = m_v / rho_glass
    # Do some geometry
    rad = sqrt(A_p/π)
    h_f = m_f/m_f0 * h_f0 
    h_d = h_f0 - h_f
    # Heat transfer from shelf
    K_shf = K_shf_f(pch(t))
    Q_shf = K_shf*A_p*(Tsh(t)-T_f) 
    Q_shw = K_shf*(A_v-A_p)*(Tsh(t)-T_f)
    # Evaluate mass flow; positive means drying is progressing. Not forced to be positive
    mflow = A_p/Rp(h_d)*(calc_psub(T_f) - pch(t)) # g/s
    Q_sub = mflow*ΔHsub # Sublimation
    # Evaluate heat transfer from wall
    Bi = uconvert(NoUnits, K_vwf*rad/k_dry)
    Q_vwf = 2π*(K_vwf*rad*h_f + k_dry*(h_f0-h_f)*S_interp(Bi)) * (T_vw-T_f)
    # Volumetric heating
    Qppp_RF_f  = 2*pi*f_RF*e_0*epp_f(T_f, f_RF)*P_per_vial(t)*B_f # W / m^3
    Qppp_RF_vw = 2*pi*f_RF*e_0*epp_vw*P_per_vial(t)*B_vw # W / m^3
    Q_RF_f = Qppp_RF_f*A_p*h_f # W
    Q_RF_vw = Qppp_RF_vw*V_vial # W
    # Check that total volumetric heating is less than input power
    if Q_RF_f + Q_RF_vw > P_per_vial(t) && t == 0u"hr"
        @warn "Energy balance of EM terms not satisfied." uconvert(u"W", Q_RF_f) uconvert(u"W", Q_RF_vw) P_per_vial(t)
    end
    # Evaluate derivatives
    # Desublimation is not allowed here: if we clamp mflow itself, then the DAE is unstable
    dm_f = min(0u"kg/s", -mflow/porosity)
    dT_f = (Q_RF_f-Q_sub+Q_vwf+Q_shf) / (m_f*cp_f) - T_f*dm_f/m_f
    dT_vw = (Q_RF_vw+Q_shw-Q_vwf) / (m_v*cp_v)
    # Strip units from derivatives; return all heat transfer terms
    return ustrip.([u"g/hr", u"K/hr", u"K/hr"], [dm_f, dT_f, dT_vw]), uconvert.(u"W", [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw, ])
end

# ```
# params = (   
#     (Rp, h_f0, c_solid, ρ_solution),
#     (K_shf_f, A_v, A_p),
#     (pch, Tsh, P_per_vial),
#     (m_f0, cp_f, m_v, cp_v, A_rad),
#     (f_RF, epp_f, epp_vw),
#     (K_vwf, B_f, B_vw, alpha),
# )
# ```

struct ParamObjRF{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, 
                T11, T12, T13, T14, T15, T16, T17, T18, T19,
                T20, T21, T22} 
    Rp::T1
    h_f0::T2
    c_solid::T3
    ρ_solution::T4
    K_shf_f::T5
    A_v::T6
    A_p::T7
    pch::T8
    Tsh::T9
    P_per_vial::T10
    m_f0::T11
    cp_f::T12
    m_v::T13
    cp_v::T14
    A_rad::T15
    f_RF::T16
    epp_f::T17
    epp_vw::T18
    K_vwf::T19
    B_f::T20
    B_vw::T21
    alpha::T22
end

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
        return (po.Rp, po.h_f0, po.c_solid, po.ρ_solution)
    elseif i == 2
        return (po.K_shf_f, po.A_v, po.A_p)
    elseif i==3 
        return (po.pch, po.Tsh, po.P_per_vial)
    elseif i == 4
        return (po.m_f0, po.cp_f, po.m_v, po.cp_v, po.A_rad)
    elseif i == 5
        return (po.f_RF, po.epp_f, po.epp_vw)
    elseif i == 6
        return (po.K_vwf, po.B_f, po.B_vw, po.alpha)
    else
        error(BoundsError, "Attempt to access LyoPronto.ParamsObjRF at index $i. Only indices 1 to 6 allowed")
    end
end