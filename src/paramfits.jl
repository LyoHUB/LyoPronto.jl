export gen_sol_conv_dim, obj_tT_conv

@doc raw"""
    gen_sol_conv_dim(KRp_prm, otherparams, u0, tspan; kwargs...)

Solve the Pikal model for primary drying with given Kv & Rp, returning the solution object and the set of parameters passed to `solve`.

- `otherparams` contains: `(hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)`
- `u0` is given as floats (not Unitful quantities), with dimensions [cm, K].
- KRp_prm has the form `[Kv, R0, A1, A2]`; this function assigns the units `[W/m^2/K, cm^2*Torr*hr/g, cm*Torr*hr/g, 1/cm]` to those numbers
- `kwargs` is passed directly (as is) to the ODE `solve` call.

This is used as a helper for [`obj_tT_conv`](@ref) to assemble a function which takes [Kv, R0, A1, A2] and returns sum squared error against experiment.
To that end, use this to make an anonymous function of the form `gen_sol = x->gen_sol_conv_dim(x, case_params, case_u0, case_tspan)`.
That function is what you will pass to `obj_tT_conv`.
"""
function gen_sol_conv_dim(KRp_prm, otherparams, u0, tspan; kwargs...)
    hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh = otherparams
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    Kshf_g = p->KRp_prm[1]*u"W/m^2/K"
    Rp_g = RpFormFit((KRp_prm[2:4].*Rp_un)...)
    new_params = ((Rp_g, hf0, c_solid, ρ_solution), (Kshf_g, Av, Ap),
                    (pch, Tsh))
    newprob = ODEProblem(LyoPronto.lyo_1d_dae_f, u0, tspan, new_params)
    sol = solve(newprob, Rodas4P(); callback=end_drying_callback, kwargs...)
    return sol, new_params
end

@doc raw"""
    obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, verbose = true)

Evaluate an objective function which compares model solution with `KRp_prm` to experimental data in `tTdat`.

Arguments:
- `gen_sol` is a function taking `[Kv, R0, A1, A2]` and returning a solution to Pikal model; see [`gen_sol_conv_dim`](@ref LyoPronto.gen_sol_conv_dim).
- `tTdat` is experimental temperature series, of the form `(time, Tf)`.
- `tweight` gives the weighting (in `K^2/hr^2`) of the end of drying in the objective, as compared to the temperature error.
- `t_end` has a default value of `0.0u"hr"`, which (if left at default) is replaced with the last time point.
- `verbose` defaults to `true`, in which case each call to this function prints info on the passed parameters, etc.
"""
function obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, verbose=true)
    if any(KRp_prm .< 0)
        return NaN
    end
    tdat, Tdat = tTdat
    if t_end == 0u"hr"
        t_end = tdat[end]
    end
    sol = gen_sol(KRp_prm)[1]
    tmd = sol.t[end] *u"hr"
    trim = tdat .< tmd
    t_trim = tdat[trim]
    if sol.u[end][1] > 1e-5
        @info "loss call: failed integration" sol
        return NaN
    end

    Tmd = sol(ustrip.(u"hr", t_trim), idxs=2).*u"K" 
    Tobj = ustrip(u"K^2", sum(abs2.(Tdat[trim] .- Tmd))/length(t_trim))
    tobj = (ustrip(u"hr", t_end[end] - tmd))^2
    verbose && @info "loss call" KRp_prm t_end-tmd extrema(Tdat[trim] .- Tmd) tobj Tobj 
    return Tobj + tweight*tobj 
end

@doc raw"""
    gen_sol_rf_dim(fitprm, params_bunch, u0, tspan; kwargs...)

Solve the lumped-capacitance model for microwave-assisted primary drying with given fit parameters, returning the solution object and the set of parameters passed to `solve`.

- `fitprm` has the form `[α, Kvwf, Bf, Bvw]`; this function assigns the units `[cm^1.5, cal/s/K/cm^2, Ω/m^2, Ω/m^2]` to those numbers.
- `params_bunch` contains a full listing of parameters used for the model, according to [`lumped_cap_rf`](@ref LyoPronto.lumped_cap_rf), including dummy values for the fit parameters.
- `u0` is given as floats (not Unitful quantities), with dimensions [g, K, K].
- `kwargs` is passed directly (as is) to the ODE `solve` call.

This is used as a helper for [`obj_tT_rf`](@ref), [`obj_tTT_rf`](@ref), [`obj_ttTT_rf`](@ref) to assemble a function which takes [α, Kvwf, Bf, Bvw] and returns sum squared error against experiment.
To that end, use this to make an anonymous function of the form `gen_sol = x->gen_sol_rf_dim(x, case_params, case_u0, case_tspan)`.
That function is what you will pass to `obj_tT_rf` (or similar).
"""
function gen_sol_rf_dim(fitprm, params_bunch, u0, tspan; kwargs...)
    prm_un = [u"cm^1.5", u"cal/s/K/cm^2", u"Ω/m^2", u"Ω/m^2"]
    newp = deepcopy(params_bunch)
    newp[end] = tuple((fitprm .* prm_un)...)
    newp = tuple(newp...)

    newprob = ODEProblem(lumped_cap_rf, u0, tspan, newp)
    sol = solve(newprob, Rodas3(autodiff=false); callback=end_drying_callback, kwargs...)
    return sol, newp
end


@doc raw"""
    obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u"hr", tweight=1, verbose=true)

Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTTdat`.

- `gen_sol` is a function taking [α, Kvwf, Bf, Bvw] and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
- `tTTdat` is experimental temperature series, of the form `(time, Tf, Tvw)`, so with frozen and vial wall temperatures taken at the same time points. 
    See also [`obj_tT_rf`](@ref) and [`obj_ttTT_rf`](@ref).
- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
- `t_end` has a default value of `0.0u"hr"`, which (if left at default) is replaced with the last given time point.
"""
function obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u"hr", tweight=1, verbose=true)
    # if fitprm[2] < 0
    if any(fitprm .< 0)
        return NaN
    end
    tloc, Tfdat, Tvwdat = tTTdat
    if t_end == 0u"hr"
        t_end = tloc[end]
    end
    sol = gen_sol(fitprm)[1]
    tmd = sol.t[end].*u"hr"
    trim = tloc .< tmd
    t_trim = tloc[trim]

    Tfmd = sol(ustrip.(u"hr", t_trim), idxs=2).*u"K"# .- 273.15
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tfmd .< 0u"K")# Ugly fix for bad interpolated values
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        @info "bad interpolation" subzero Tfmd[subzero]
    end
    Tvwmd = sol(ustrip.(u"hr", t_trim), idxs=3).*u"K"# .- 273.15
    Tfobj = sum(abs2.(Tfdat[trim] .- Tfmd))/length(t_trim)
    Tvwobj = sum(abs2.(Tvwdat[trim] .- Tvwmd))/length(t_trim)
    tobj = ((t_end - tmd))^2
    verbose && @info "loss call" fitprm tmd t_end Tfobj Tvwobj 
    return Tfobj/u"K^2" + Tvwobj/u"K^2" + tweight*tobj/u"hr^2"
end


@doc raw"""
    obj_tT_rf(fitprm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, Tvw_end = 0.0u"K", verbose=true)

Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTdat`.

- `gen_sol` is a function taking `[α, Kvwf, Bf, Bvw]` and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
- `tTTdat` is experimental temperature series, of the form `(time, Tf)`, so with frozen temperatures only. 
    See also [`obj_tTT_rf`](@ref) and [`obj_ttTT_rf`](@ref).
- `Tvw_end` defaults to `0.0u"K"`, in which case vial wall temperatures are excluded from the objective. 
    If another value is passed, then the final model vial wall temperature is compared to that value and included in the objective.
- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
- `t_end` has a default value of `0.0u"hr"`, which (if left at default) is replaced with the last given time point.
"""
function obj_tT_rf(fitprm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, Tvw_end = 0.0u"K", verbose=true)
    if any(fitprm .< 0)
        return NaN
    end
    tloc, Tfdat = tTdat
    if t_end == 0u"hr"
        t_end = tloc[end]
    end
    sol = gen_sol(fitprm)[1]
    tmd = sol.t[end].*u"hr"
    trim = tloc .< tmd
    t_trim = tloc[trim]

    Tfmd = sol(ustrip.(u"hr", t_trim), idxs=2).*u"K"
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tfmd .< 0u"K")  
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        @info "bad interpolation replaced" subzero Tfmd[subzero]
    end
    Tfobj = sum(abs2.(Tfdat[trim] .- Tfmd))/length(t_trim)
    if Tvw_end > 0u"K"
        Tvw_obj = (sol[3, end]*u"K" - Tvw_end)^2
    else
        Tvw_obj = 0u"K^2"
    end
    tobj = ((t_end - tmd))^2
    verbose && @info "loss call" fitprm tmd t_end Tfobj Tvw_obj
    return Tfobj/u"K^2" + Tvw_obj/u"K^2" + tweight*tobj/u"hr^2"
end

@doc raw"""
    obj_ttTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u"hr", tweight=1, verbose=true)

Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTdat`.

- `gen_sol` is a function taking `[α, Kvwf, Bf, Bvw]` and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
- `tTTdat` is experimental temperature series, of the form `(time_Tf, time_vw, Tf, Tvw)`, so with Tf and Tvw having separate time points. 
    This is useful if there is an early temperature rise in Tf, but Tvw continues to be reliable, so the model can fit to as much of Tvw as reasonable.
    See also [`obj_tTT_rf`](@ref) and [`obj_tT_rf`](@ref).
- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
- `t_end` has a default value of `0.0u"hr"`, which (if left at default) is replaced with the last given time point.
"""
function obj_ttTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u"hr", tweight=1, verbose=true)
    if any(fitprm .< 0)
        return NaN
    end
    tf, tvw, Tfdat, Tvwdat = tTTdat
    if t_end == 0u"hr"
        t_end = max(tf[end], tvw[end])
    end
    sol = gen_sol(fitprm)[1]
    tmd = sol.t[end].*u"hr"
    ftrim = tf .< tmd
    vwtrim = tvw .< tmd
    tf_trim = tf[ftrim]
    tvw_trim = tvw[vwtrim]

    Tfmd = sol(ustrip.(u"hr", tf_trim), idxs=2).*u"K"
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tfmd .< 0u"K")
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        @info "bad interpolation" subzero Tfmd[subzero]
    end
    Tvwmd = sol(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
    Tfobj = sum(abs2.(Tfdat[ftrim] .- Tfmd))/length(tf_trim)
    Tvwobj = sum(abs2.(Tvwdat[vwtrim] .- Tvwmd))/length(tvw_trim)
    tobj = ((t_end - tmd))^2
    verbose && @info "loss call" fitprm tmd t_end Tfobj Tvwobj 
    return Tfobj/u"K^2" + Tvwobj/u"K^2" + tweight*tobj/u"hr^2"
end
