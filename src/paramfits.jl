export gen_sol_conv_dim, obj_tT_conv

"""
    obj_tT_conv(KRp_prm, otherparams, u0, tTdat; tweight=1)

Evaluate an objective function which compares model solution with `KRp_prm` to experimental data in `tTdat`.

- `otherparams` contains: `(hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)`
- `u0` is unitless, but with dimensions [cm, K].
- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
"""
function gen_sol_conv_dim(KRp_prm, otherparams, u0, tspan; kwargs...)
    hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh = otherparams
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    Kshf_g = p->KRp_prm[1]*u"W/m^2/K"
    Rp_g = RpFormFit((KRp_prm[2:4].*Rp_un)...)
    new_params = ((Rp_g, hf0, c_solid, ρ_solution), (Kshf_g, Av, Ap),
                    (pch, Tsh))
    newprob = ODEProblem(LyoPronto.lyo_1d_dae_f, u0, tspan, new_params)
    # @info "getting sol" Rp
    sol = solve(newprob, Rodas4P(); callback=end_drying_callback, kwargs...)
    return sol, new_params
end

"""
    obj_tT_conv(KRp_prm, gen_sol, tTdat; tweight=1)

Evaluate an objective function which compares model solution with `KRp_prm` to experimental data in `tTdat`.

- `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
"""
function obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1)
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
    @info "loss call" KRp_prm t_end-tmd extrema(Tdat[trim] .- Tmd) tobj Tobj 
    return Tobj + tweight*tobj 
end

function gen_sol_rf_dim(fitprm, params_bunch, u0, tspan; kwargs...)
    prm_un = [u"cm^1.5", u"cal/s/K/cm^2", u"Ω/m^2", u"Ω/m^2"]
    newp = deepcopy(params_bunch)
    newp[end] = tuple((fitprm .* prm_un)...)
    newp = tuple(newp...)

    newprob = ODEProblem(lumped_cap_rf, u0, tspan, newp)
    sol = solve(newprob, Rodas3(autodiff=false); callback=end_drying_callback, kwargs...)
    return sol, newp
end


function obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=0u"hr", tweight=1)
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
    if any(Tfmd .< 0u"K")# Ugly fix for bad interpolated values
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        @info "bad interpolation" subzero Tfmd[subzero]
    end
    Tvwmd = sol(ustrip.(u"hr", t_trim), idxs=3).*u"K"# .- 273.15
    Tfobj = sum(abs2.(Tfdat[trim] .- Tfmd))/length(t_trim)
    Tvwobj = sum(abs2.(Tvwdat[trim] .- Tvwmd))/length(t_trim)
    tobj = ((t_end - tmd))^2
    @info "loss call" fitprm tmd t_end Tfobj Tvwobj 
    return Tfobj/u"K^2" + Tvwobj/u"K^2" + tweight*tobj/u"hr^2"
end


function obj_tT_rf(fitprm, gen_sol, tTdat; t_end=0u"hr", tweight=1, Tvw_end = 0.0u"K")
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

    Tfmd = sol(ustrip.(u"hr", t_trim), idxs=2).*u"K"# .- 273.15
    if any(Tfmd .< 0u"K")# Ugly fix for bad interpolated values
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        @info "bad interpolation" subzero Tfmd[subzero]
    end
    Tfobj = sum(abs2.(Tfdat[trim] .- Tfmd))/length(t_trim)
    if Tvw_end > 0u"K"
        Tvw_obj = (sol[3, end]*u"K" - Tvw_end)^2
    else
        Tvw_obj = 0u"K^2"
    end
    tobj = ((t_end - tmd))^2
    @info "loss call" fitprm tmd t_end Tfobj Tvw_obj
    return Tfobj/u"K^2" + Tvw_obj/u"K^2" + tweight*tobj/u"hr^2"
end

function obj_ttTT_rf(fitprm, gen_sol, tTTdat; t_end=0u"hr", tweight=1)
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

    Tfmd = sol(ustrip.(u"hr", tf_trim), idxs=2).*u"K"# .- 273.15
    if any(Tfmd .< 0u"K")# Ugly fix for bad interpolated values
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        @info "bad interpolation" subzero Tfmd[subzero]
    end
    Tvwmd = sol(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
    Tfobj = sum(abs2.(Tfdat[ftrim] .- Tfmd))/length(tf_trim)
    Tvwobj = sum(abs2.(Tvwdat[vwtrim] .- Tvwmd))/length(tvw_trim)
    tobj = ((t_end - tmd))^2
    @info "loss call" fitprm tmd t_end Tfobj Tvwobj 
    return Tfobj/u"K^2" + Tvwobj/u"K^2" + tweight*tobj/u"hr^2"
end
