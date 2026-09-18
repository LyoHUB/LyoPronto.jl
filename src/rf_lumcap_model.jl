
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
# TODO: this line takes 1.5s to precompile, consider computing elsewhere 
# If done when params are constructed, will incur the cost repeatedly in fitting, which is worse.
const S_samp = shapefac.(Bi_samp)
const S_interp = LinearInterpolation(S_samp, Bi_samp, extrapolation=ExtrapolationType.Linear)

const RF_PARAMS_DOC = """

```
params = ParamObjRF((   
    (Rp, hf0, cSolid, ρsolution),
    (Kshf_f, Av, Ap),
    (pch, Tsh, P_per_vial),
    (mf0, cpf, mv, cpv),
    (f_RF, eppf, eppvw),
    (Kvwf, Bf, Bvw),
))
```

The parameters should all be Unitful quantities with appropriate dimensions, with some exceptions which are callables returning quantities.
See [`RpFormFit`](@ref) and [`RampedVariable`](@ref) for convenience types that can help with these cases.
- `Rp(x)` with `x` a length returns mass transfer resistance (as a Unitful quantity)
- `Kshf_f(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)`, `P_per_vial(t)` return shelf temperature, chamber pressure, and microwave power respectively at time `t`.

- `Arad` and `alpha` were used in a prior version of the model, and are not used in the 
    current version; they will be removed in a future version.
"""


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
    f_RF
    eppf
    eppvw
    Kvwf
    Bf
    Bvw
end

@doc """
    $(TYPEDEF)

The `ParamObjRF` type is a container for the parameters used in the RF model.


Since it has many fields, the recommended constructor accepts a tuple of tuples, 
as follows, to help avoid ordering mistakes:

$(RF_PARAMS_DOC)
"""
ParamObjRF

function ParamObjRF(tuptup::Tuple) 
    # Check the proper number of fields
    length.(tuptup) == (4, 3, 3, 4, 3, 3) || error("ParamObjRF tuple-of-tuple structure is wrong.")
    # Create the object
    po = ParamObjRF(tuptup[1]..., tuptup[2]...,
                tuptup[3]..., tuptup[4]...,
                tuptup[5]..., tuptup[6]...,)
    # Check that callables return proper outputs
    po.Kshf(1u"Torr") * u"m^2"*u"K" isa Unitful.Power || error("Kshf does not return heat transfer coeff")
    # Note that, odd though it sounds, Rp does have dimensions of velocity
    po.Rp(1u"cm") isa Unitful.Velocity || error("Rp does not return a mass transfer resistance")
    po.pch(1.0u"hr") isa Unitful.Pressure || error("pch does not return a pressure")
    po.Tsh(1.0u"hr") isa Unitful.Temperature || error("Tsh does not return an absolute temperature")
    return po 
end


"""
    $(SIGNATURES)

Compute the mass flow and  heat transfer terms for the lumped-capacitance microwave-assisted model.

Returns a named tuple with the following fields, all as Unitful quantities:
- `md`: mass flow rate (g/hr)
- `Q_shf`: heat transfer from shelf to product (W)
- `Q_vwf`: heat transfer from vial wall to product (W)
- `Q_RF_f`: volumetric heating of product (W)
- `Q_RF_vw`: volumetric heating of vial wall (W)
- `Q_shw`: heat transfer from shelf to vial wall (W)
"""
@inline function calc_md_Q(u, po::ParamObjRF, tn)
    # Unpack all the parameters
    (;Rp, hf0, csolid, ρsolution,
    Kshf, Av, Ap,
    pch, Tsh, P_per_vial, 
    mf0, mv,
    f_RF, eppf, eppvw,
    Kvwf, Bf, Bvw) = po
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
    Q_shf = Kshft*Ap*(Tsh(t)-T_f) |> u"W"
    Q_shw = Kshft*(Av-Ap)*(Tsh(t)-T_vw) |> u"W"
    # Evaluate mass flow; positive means drying is progressing. Not forced to be positive
    mflow = Ap/Rp(h_d)*(calc_psub(T_f) - pch(t)) # g/s
    # Evaluate heat transfer from wall
    # TODO: consider precalculating Bi and shape factor in ParamObjRF constructor
    Bi = uconvert(NoUnits, Kvwf*rad/k_dry)
    Q_vwf = 2π*(Kvwf*rad*h_f + k_dry*(hf0-h_f)*S_interp(Bi)) * (T_vw-T_f) |> u"W"
    # Volumetric heating
    Qppp_RF_f  = 2*pi*f_RF*e_0*eppf(T_f, f_RF)*P_per_vial(t)*Bf # W / m^3
    Qppp_RF_vw = 2*pi*f_RF*e_0*eppvw*P_per_vial(t)*Bvw # W / m^3
    Q_RF_f = Qppp_RF_f*Ap*h_f |> u"W" # W
    Q_RF_vw = Qppp_RF_vw*V_vial |> u"W"# W
    Q_sub = mflow*ΔHsub # Sublimation
    # Check that total volumetric heating is less than input power
    if Q_RF_f + Q_RF_vw > P_per_vial(t) && t == 0u"hr"
        @warn "Energy balance of EM terms not satisfied." Q_RF_f Q_RF_vw P_per_vial(t)
    end
    return (; md=mflow, Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw)
end

"""
    $(SIGNATURES)

Compute the right-hand-side function for the ODEs making up the lumped-capacitance microwave-assisted model.

To access the values of the various heat transfer terms, use `[calc_md_Q](@ref)` to compute them; that function is used internally by this function.

`du` refers to `[dmf/dt, dTf/dt, dTvw/dt]`, with `u = [mf, Tf, Tvw]`.
`u` is taken without units but assumed to have the units of `[g, K, K]` (which is internally added).
`tn` is assumed to be in hours (internally added), so `dudt` is returned with assumed units `[g/hr, K/hr, K/hr]` to be consistent.

Use the `ParamObjRF` type to hold the parameters. 
$(RF_PARAMS_DOC)

"""
function lumped_cap_rf!(du, u, params::ParamObjRF, tn)

    # Compute heat transfer rates
    (; md, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw) = calc_md_Q(u, params, tn)
    mflow = md

    (; csolid, ρsolution,
    cpf, mv, cpv ) = params
    # Dimensionalize the state variables
    m_f = u[1]*u"g"
    T_f = u[2]*u"K"
    porosity = (ρsolution - csolid)/ρsolution

    Q_sub = mflow*ΔHsub # Sublimation

    # Block desublimation
    dm_f = min(0.0u"kg/s", -mflow/porosity)
    dT_f =  (Q_shf+Q_vwf+Q_RF_f -Q_sub) / (m_f*cpf) - T_f*dm_f/m_f
    dT_vw = (Q_shw-Q_vwf+Q_RF_vw) / (mv*cpv)

    # Strip units from derivatives
    du[1] = ustrip(u"g/hr", dm_f)
    du[2] = ustrip(u"K/hr", dT_f)
    du[3] = ustrip(u"K/hr", dT_vw)
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
        tstops = tstops, callback=end_drying_callback, dt=0.1)
end


@doc raw"""
    qrf_integrate(sol, RF_params)

Compute the integral over time of each heat transfer mode in the lumped capacitance model.
RF_params should represent the same parameters used to generate the solution `sol`,
which (if OrdinaryDiffEq doesn't change) can likely be accessed as `sol.prob.p`.

Returns a Dict{String, Quantity{...}}, with string keys `Qsub, Qshf, Qvwf, QRFf, QRFvw, Qshw`.
"""
function qrf_integrate(sol, RF_params::ParamObjRF)

    # Using an IntegratingSumCallback would be more elegant, but at last attempt
    # it struggled with unitful values in the arrays.
    # So we do a manual Riemann integration on the solution output
    history = Table(map(sol.t) do ti
        calc_md_Q(sol(ti), RF_params, ti)
    end)

    t = sol.t*u"hr"
    weights = fill(first(t), length(sol.t))
    dt = diff(t)
    weights[begin:end-1] += dt./2
    weights[begin+1:end] += dt./2

    names = [:Q_sub, :Q_shf, :Q_vwf, :Q_RF_f, :Q_RF_vw, :Q_shw]
    qinteg = map(names) do q
        sum(getproperty.(history, q) .* weights) |> u"W*hr"
    end
    # TODO: consider returning differently
    return Dict("Qsub"=>qinteg[1], 
                "Qshf"=>qinteg[2],
                "Qvwf"=>qinteg[3],
                "QRFf"=>qinteg[4],
                "QRFvw"=>qinteg[5],
                "Qshw"=>qinteg[6])
end