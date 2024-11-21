export gen_sol_conv_dim, obj_tT_conv
export obj_expT, genobj_posprm

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
function gen_sol_conv_dim(KRp_prm, po::ParamObjPikal, u0, tspan; kwargs...)
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    Kshf_g = p->KRp_prm[1]*u"W/m^2/K"
    Rp_g = RpFormFit((KRp_prm[2:4].*Rp_un)...)
    new_params = @set po.Rp = Rp_g
    @reset new_params.Kshf = Kshf_g
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
- `t_end` defaults to `missing`, in which case it is excluded from the objective.
    If another value is passed, then the model drying time is compared to that value and included in the objective.
- `verbose` defaults to `true`, in which case each call to this function prints info on the passed parameters, etc.
"""
function obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=missing, tweight=1, verbose=true)
    if any(KRp_prm .< 0)
        return NaN
    end
    tdat, Tdat = tTdat
    sol = gen_sol(KRp_prm)[1]
    if sol.u[end][1] > 1e-5
        @info "loss call: failed integration" sol
        return NaN
    end
    pdfit = PrimaryDryFit(tdat, (Tdat,), t_end)
    verbose && @info "gen_sol" KRp_prm
    return obj_expT(sol, pdfit; tweight=tweight, verbose=true)
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
    prm_un = [u"cal/s/K/cm^2", u"Ω/m^2", u"Ω/m^2", u"cm^1.5"]
    newp = deepcopy(params_bunch)
    newp[end] = tuple((fitprm .* prm_un)...)
    newp = tuple(newp...)
    newprob = ODEProblem(lumped_cap_rf, u0, tspan, newp)
    sol = solve(newprob, Rodas3(autodiff=false); callback=end_drying_callback, kwargs...)
    return sol, newp
end
function gen_sol_rf_dim(fitprm, po::ParamObjRF, u0, tspan; kwargs...)
    prm_un = [u"cal/s/K/cm^2", u"Ω/m^2", u"Ω/m^2", u"cm^1.5"]
    newp = @set po.alpha = fitprm[4]*prm_un[4]
    @reset newp.K_vwf = fitprm[1]*prm_un[1]
    @reset newp.B_f = fitprm[2]*prm_un[2]
    @reset newp.B_vw = fitprm[3]*prm_un[3]
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
- `t_end` has a default value of `0.0u"hr"`, which (if left at default) is not included in the objective function.
- `t_end` defaults to `missing`, in which case it is excluded from the objective.
    If another value is passed, then the model drying time is compared to that value and included in the objective.
"""
function obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=missing, tweight=1, verbose=true)
    if any(fitprm .< 0)
        return NaN
    end
    tloc, Tfdat, Tvwdat = tTTdat
    pdfit = PrimaryDryFit(tloc, (Tfdat,), tloc, (Tvwdat,), t_end)
    sol = gen_sol(fitprm)[1]
    verbose && @info "gen_sol" fitprm 
    return obj_expT(sol, pdfit; tweight=tweight, verbose=verbose)
end


@doc raw"""
    obj_tT_rf(fitprm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, Tvw_end = 0.0u"K", verbose=true)

Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTdat`.

- `gen_sol` is a function taking `[α, Kvwf, Bf, Bvw]` and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
- `tTTdat` is experimental temperature series, of the form `(time, Tf)`, so with frozen temperatures only. 
    See also [`obj_tTT_rf`](@ref) and [`obj_ttTT_rf`](@ref).
- `Tvw_end` defaults to `missing`, in which case vial wall temperatures are excluded from the objective. 
    If another value is passed, then the final model vial wall temperature is compared to that value and included in the objective.
- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
- `t_end` defaults to `missing`, in which case it is excluded from the objective.
    If another value is passed, then the model drying time is compared to that value and included in the objective.
"""
function obj_tT_rf(fitprm, gen_sol, tTdat; t_end=missing, tweight=1, Tvw_end = missing, verbose=true)
    if any(fitprm .< 0)
        return NaN
    end
    tloc, Tfdat = tTdat
    pdfit = PrimaryDryFit(tloc, (Tfdat,), missing, Tvw_end, t_end)
    sol = gen_sol(fitprm)[1]
    verbose && @info "gen_sol" fitprm 
    return obj_expT(sol, pdfit; tweight=tweight, verbose=verbose)
end

@doc raw"""
    obj_ttTT_rf(fitprm, gen_sol, ttTTdat; t_end=missing, tweight=1, verbose=true)

Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTdat`.

- `gen_sol` is a function taking `[α, Kvwf, Bf, Bvw]` and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
- `ttTTdat` is experimental temperature series, of the form `(time_Tf, time_vw, Tf, Tvw)`, so with Tf and Tvw having separate time points. 
    This is useful if there is an early temperature rise in Tf, but Tvw continues to be reliable, so the model can fit to as much of Tvw as reasonable.
    See also [`obj_tTT_rf`](@ref) and [`obj_tT_rf`](@ref).
- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
- `t_end` defaults to `missing`, in which case it is excluded from the objective.
    If another value is passed, then the model drying time is compared to that value and included in the objective.
"""
function obj_ttTT_rf(fitprm, gen_sol, ttTTdat; t_end=missing, tweight=1, verbose=true)
    if any(fitprm .< 0)
        return NaN
    end
    tf, tvw, Tfdat, Tvwdat = ttTTdat
    pdfit = PrimaryDryingFit(tf, (Tfdat,), tvw, (Tvwdat), t_end)
    sol = gen_sol(fitprm)[1]
    verbose && @info "gen_sol" fitprm 
    return obj_expT(sol, pdfit; tweight=tweight, verbose=verbose)
end

@doc raw"""
    obj_expT(sol, pdfit; tweight=1, verbose=true, rf = true)

Experimental (in the software engineering sense)!
Evaluate an objective function which compares model solution computed by `sol` to experimental data in `pdfit`.

- `sol` is a solution to an appropriate model; see [`gen_sol_conv_dim`](@ref) and [`gen_sol_rf_dim`](@ref) for some helper functions for this.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable returned by `gen_sol` is assumed to be temperature, as is true for [`gen_sol_rf_dim`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.

I've considered writing several methods and dispatching on `pdfit` somehow, which would be cool and might individually be easier to read. But control flow might be harder to document and explain, and this should work just fine.
"""
function obj_expT(sol, pdfit; tweight=1.0, verbose = false)
    if sol.retcode !== ReturnCode.Terminated
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
        return NaN
    end
    tmd = sol.t[end].*u"hr"
    ftrim = pdfit.t_Tf .< tmd
    tf_trim = pdfit.t_Tf[ftrim]
    Tfmd = sol(ustrip.(u"hr", tf_trim), idxs=2).*u"K"
    # Compute temperature objective for all frozen temperatures
    Tfobj = mapreduce(+, pdfit.Tfs, pdfit.Tf_iend) do Tf, itf
        trim = 1:min(itf, length(tf_trim))
        return sum(abs2.(Tf[trim] .- Tfmd[trim]))/length(trim)
    end
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tfmd .< 0u"K")
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        verbose && @info "bad interpolation" subzero Tfmd[subzero]
    end
    if ismissing(pdfit.Tvws) # No vial wall temperatures
        Tvwobj = 0u"K^2"
    elseif ismissing(pdfit.t_Tvw) # Provide only an endpoint temperature
        Tvwend = pdfit.Tvws
        Tvwobj = (sol[3, end]*u"K" - uconvert(u"K", Tvwend))^2
    else # Regular case of fitting to at least one full temperature series
        vwtrim = pdfit.t_Tvw .< tmd
        tvw_trim = pdfit.t_Tvw[vwtrim]
        Tvwmd = sol(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
        # Compute temperature objective for all vial wall temperatures
        Tvwobj = mapreduce(+, pdfit.Tvws, pdfit.Tvw_iend) do Tvw, itvw
            trim = 1:min(itvw, length(tvw_trim))
            return sum(abs2.(Tvw[trim] .- Tvwmd[trim]))/length(trim)
        end
    end
    if ismissing(pdfit.t_end)
        tobj = 0u"hr^2"
    else # Compare drying time
        tobj = ((pdfit.t_end - tmd))^2
    end
    verbose && @info "loss call" tmd tobj Tfobj Tvwobj 
    return Tfobj/u"K^2" + Tvwobj/u"K^2" + tweight*tobj/u"hr^2"
end

function genobj_posprm(gen_sol, obj, fitdat; kwargs...)
    return (x->(any(x .< 0) && return NaN; obj(gen_sol(x)[1], fitdat; kwargs...)))
end
