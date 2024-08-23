export lyo_1d!, lyo_1d_dae_f, subflux_Tsub

const ΔHsub = 678u"cal/g"
const k_ice = 2.45u"W/m/K"

end_cond(u, t, integ) = u[1] - 1e-10 # When reaches 1e-10, is basically zero
const end_drying_callback = ContinuousCallback(end_cond, terminate!)

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

lyo_1d_dae_f = ODEFunction(lyo_1d_dae!, mass_matrix=lyo_1d_mm)

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