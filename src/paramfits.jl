export gen_sol_KRp, gen_sol_Rp, gen_sol_rf
export obj_KRp, obj_Rp, obj_KBB
export obj_expT, err_expT 

export copy_nt_into_struct
function copy_nt_into_struct(po, nt::NamedTuple)
    @warn "Copying arbitrary NamedTuple into struct. Type unstable. If doing this repeatedly, define a new method for copy_nt_into_struct" nt
    po_new = deepcopy(po)
    # Merge the parameters in the NamedTuple into the ParamObjPikal object
    for (k, v) in zip(keys(nt), values(nt))
    # map(keys(nt), values(nt)) do k, v
        # @info "copying" k v PropertyLens(k)
        po_new = set(po_new, PropertyLens{k}(), v)
        nothing
    end
    return po_new
end
function copy_nt_into_struct(po, nt::NamedTuple{(:R0, :A1, :A2)})
    po_new = deepcopy(po)
    @reset po_new.Rp.R0 = nt.R0
    @reset po_new.Rp.A1 = nt.A1
    @reset po_new.Rp.A2 = nt.A2
end
function copy_nt_into_struct(rp::RpFormFit, nt::NamedTuple{(:R0, :A1, :A2)})
    rp_new = deepcopy(rp)
    @reset rp_new.R0 = nt.R0
    @reset rp_new.A1 = nt.A1
    @reset rp_new.A2 = nt.A2
end
function copy_nt_into_struct(po, nt::NamedTuple{(:Kshf, :R0, :A1, :A2)})
    po_new = deepcopy(po)
    @reset po_new.Kshf = nt.Kshf
    @reset po_new.Rp.R0 = nt.R0
    @reset po_new.Rp.A1 = nt.A1
    @reset po_new.Rp.A2 = nt.A2
end

"""
    $(SIGNATURES)

Evaluate an objective function which compares model solution computed by `sol` to experimental data in `pdfit`.

- `sol` is a solution to an appropriate model; see [`gen_sol_Rp`](@ref), [`gen_sol_KRp`](@ref), and [`gen_sol for some helper functions for this.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for [`gen_sol_rf`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.

I've considered writing several methods and dispatching on `pdfit` somehow, which would be cool and might individually be easier to read. But control flow might be harder to document and explain, and this should work just fine.
"""
function obj_expT(sol::ODESolution, pdfit::PrimaryDryFit{TT1, TT2, TT3, TT4, TT5, TTvw, TTvwi, Tte}, ; tweight=1.0, verbose = false) where {TT1, TT2, TT3, TT4, TT5, TTvw, TTvwi, Tte}
    if sol.retcode !== ReturnCode.Terminated
        # hf_end = sol[1, end]*u"cm"
        verbose && @warn "ODE solve did not reach end of drying. Either parameters are bad, or tspan is not large enough." sol.retcode sol.prob.p.hf0 sol[end]
        return NaN
    end
    tmd = sol.t[end].*u"hr"
    nt = length(sol.t) - 1
    preinterp = all(sol.t[begin:nt÷2]*u"hr" .≈ pdfit.t[begin:nt÷2]) # If ODE solution has been interpolated already to time points...
    # Compute temperature objective for all frozen temperatures
    if preinterp
        Tfmd = sol[2, begin:end-1].*u"K" # Leave off last time point because is end time
    else
        ftrim = pdfit.t .< tmd
        tf_trim = pdfit.t[ftrim]
        Tfmd = sol(ustrip.(u"hr", tf_trim), idxs=2).*u"K"
        # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
        # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
        if any(Tfmd .< 0u"K")
            subzero = findall(Vector(Tfmd .< 0u"K"))
            Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
            verbose && @info "bad interpolation" subzero Tfmd[subzero]
        end
    end
    Tfobj = 0.0u"K^2"
    for (j, iend) in enumerate(pdfit.Tf_iend)
        trim = min(iend, length(Tfmd))
        Tfobj += sum(abs2, (pdfit.Tfs[j][begin:trim] .- Tfmd[begin:trim]))/trim
    end
    if TTvw == Missing # No vial wall temperatures, encoded in type
        Tvwobj = 0u"K^2"
    elseif TTvwi == Missing # Endpoint only temperature, encoded in type
        Tvwend = pdfit.Tvws
        Tvwobj = (sol[3, end]*u"K" - uconvert(u"K", Tvwend))^2
    else # Regular case of fitting to at least one full temperature series
        if preinterp
            Tvwmd = sol[3, begin:end-1].*u"K" # Leave off last time point because is end time
        else
            vwtrim = pdfit.t .< tmd
            tvw_trim = pdfit.t[vwtrim]
            Tvwmd = sol(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
        end
        # Compute temperature objective for all vial wall temperatures
        Tvwobj = 0.0u"K^2"
        for (j, iend) in enumerate(pdfit.Tvw_iend)
            trim = min(iend, length(Tvwmd))
            Tvwobj += sum(abs2, (pdfit.Tvws[j][begin:trim] .- Tvwmd[begin:trim]))/trim
        end
    end
    if Tte == Missing # No drying time provided: encoded in type
        tobj = 0.0u"hr^2"
    else # Compare drying time
        tobj = ((pdfit.t_end - tmd))^2
    end
    verbose && @info "loss call" tmd tobj Tfobj Tvwobj 
    return ustrip(u"K^2", Tfobj + Tvwobj) + tweight*ustrip(u"hr^2", tobj)
end

function obj_expT(sol, pdfit; verbose=false, kwargs...)
    verbose && @warn "`obj_expT` got passed improper args. Might not be a problem, but check." sol
    if isnan(sol)
        return NaN
    end
    error("Improper call to `obj_expT`.")
end

@doc raw"""
    err_expT(sol, pdfit; tweight=1, verbose=true, rf = true)

Evaluate the error between model solution `sol` to experimental data in `pdfit`.

In contrast to `obj_expT()`, this function returns an array of all the errors, which would be squared and summed to produce an objective function.

- `sol` is a solution to an appropriate model; see [`gen_sol_Rp`](@ref), [`gen_sol_KRp`](@ref), and [`gen_sol_rf`](@ref) for some helper functions for this.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
Each time series, plus the end time, is given equal weight by dividing by its length; error is given in K (but `ustrip`ped).

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for [`gen_sol_rf`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.
"""
function err_expT(sol, pdfit; tweight=1.0, verbose = false)
    if sol.retcode !== ReturnCode.Terminated
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
        return NaN
    end
    tmd = sol.t[end].*u"hr"
    ftrim = pdfit.t .< tmd
    tf_trim = pdfit.t[ftrim]
    Tfmd = sol(ustrip.(u"hr", tf_trim), idxs=2).*u"K"
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tfmd .< 0u"K")
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        verbose && @info "bad interpolation" subzero Tfmd[subzero]
    end
    # Initialize error array with frozen temperatures
    # Compute temperature errors for all frozen temperatures
    errs = mapreduce(vcat, pdfit.Tfs, pdfit.Tf_iend) do Tf, itf
        trim = 1:min(itf, length(tf_trim))
        errs = ustrip.(u"K", Tf[trim] .- Tfmd[trim]) ./ sqrt(itf) # Weight errors by number of data points here
        extra = fill(0.0, itf - length(errs))
        vcat(errs, extra)
    end

    # Concatenate end time to array, if present
    if !ismissing(pdfit.t_end)
        t_err = (pdfit.t_end - tmd)
        push!(errs, ustrip(u"hr", t_err))
    end
    # If present, vcat vial wall temperatures
    if !ismissing(pdfit.Tvws) # At least one vial wall temperature
        if ismissing(pdfit.Tvw_iend) # Only an endpoint temperature provided
            Tvwend = pdfit.Tvws
            Tvw_err = sol[3, end]*u"K" - uconvert(u"K", Tvwend)
            push!(errs, ustrip(u"K", Tvw))
        else # Regular case of fitting to at least one full temperature series
            vwtrim = pdfit.t .< tmd
            tvw_trim = pdfit.t[vwtrim]
            Tvwmd = sol(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
            # Compute temperature objective for all vial wall temperatures
            Tvwobj = mapreduce(+, pdfit.Tvws, pdfit.Tvw_iend) do Tvw, itvw
                trim = 1:min(itvw, length(tvw_trim))
                Tvw_err = Tvw[trim] .- Tvwmd[trim]
                push!(errs, ustrip(u"K", Tvw_err / sqrt(itvw)))
                push!(errs, fill(0.0, itvw - length(Tvw_err)))
            end
        end
    end
    verbose && @info "loss call" tmd size(errs) sum(abs2.(errs))
    return errs
end

@doc raw"""
    gen_sol_conv(fitlog, tr, po::ParamObjPikal[, fitdat]; saveat=[], kwargs...)

Simulate primary drying with the Pikal model, given a vector of guesses, a mapping `tr` from `fitlog` to named coefficients, and other parameters in `po`.

`tr` should be a `TransformTuple` object, from TransformVariables, which maps e.g. a vector of 3
real numbers to 3 a NamedTuple with `R0, A1, A2` as keys and appropriate Unitful dimensions on the values.
Inside, this function runs 
```
fitprm = transform(tr, fitlog)
new_params = copy_nt_into_struct(po, fitprm)
```
which is done here to avoid code duplication.

So, to choose which parameters to fitting, all that is necessary is to provide an appropriate transform `tr`
and add a method of `copy_nt_into_struct` for the desired parameters.
Therefore this function can be used for both K-Rp fitting, or just Rp, or just a subset of the 3 Rp coefficients.

If given, `fitdat` is used to set `saveat` for the ODE solution.

Other `kwargs` are passed directly (as is) to the ODE `solve` call.
"""
function gen_sol_pd(fitlog, tr, po; saveat=[], kwargs...)
    fitprm = transform(tr, fitlog)
    new_params = setproperties(po, fitprm)
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3(); saveat, kwargs...)
    return sol
end
function gen_sol_pd(fitlog, tr, po, fitdat; kwargs...)
    sol = gen_sol_Rp(fitlog, tr, po; saveat=ustrip.(u"hr", fitdat.t), kwargs...)
    return sol
end
"""
    obj_conv(fitlog, tpf; tweight=1.0, verbose=false)

Calculate the sum of squared error (objective function) for fitting parameters to primary drying data.
This directly calls [`gen_sol_conv`](@ref), then [`obj_expT`](@ref), so see those docstrings.

"""
function obj_pd(fitlog, tpf; tweight=1.0, verbose=false)
    rtype = eltype(Rp_log)
    sol = gen_sol_Rp(Rp_log, tpf...)
    return rtype(obj_expT(sol, tpf[3], tweight=tweight, verbose=verbose))
end



# @doc raw"""
#     obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, verbose = true)

# Evaluate an objective function which compares model solution with `KRp_prm` to experimental data in `tTdat`.

# Arguments:
# - `gen_sol` is a function taking `[Kv, R0, A1, A2]` and returning a solution to Pikal model; see [`gen_sol_conv_dim`](@ref LyoPronto.).
# - `tTdat` is experimental temperature series, of the form `(time, Tf)`.
# - `tweight` gives the weighting (in `K^2/hr^2`) of the end of drying in the objective, as compared to the temperature error.
# - `t_end` defaults to `missing`, in which case it is excluded from the objective.
#     If another value is passed, then the model drying time is compared to that value and included in the objective.
# - `verbose` defaults to `true`, in which case each call to this function prints info on the passed parameters, etc.
# """
# function obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=missing, tweight=1, verbose=true)
#     if any(KRp_prm .< 0)
#         return NaN
#     end
#     tdat, Tdat = tTdat
#     sol = gen_sol(KRp_prm)[1]
#     if sol.u[end][1] > 1e-5
#         @info "loss call: failed integration" sol
#         return NaN
#     end
#     pdfit = PrimaryDryFit(tdat, (Tdat,), t_end)
#     verbose && @info "gen_sol" KRp_prm
#     return obj_expT(sol, pdfit; tweight=tweight, verbose=true)
# end

const rfprm_base_scale = (1.0u"W/m^2/K", 1e7u"Ω/m^2", 1e7u"Ω/m^2")

@doc raw"""
    gen_sol_rf(KBB_log, po::ParamObjRF[, fitdat::PrimaryDryFit]; rfprm_scale = rfprm_base_scale, kwargs...)

Solve the lumped-capacitance model for microwave-assisted primary drying with given fit parameters, returning the solution object and the set of parameters passed to `solve`.

- `fitprm` has the form `log.([Kvwf, Bf, Bvw])`; this function takes the exponential, then multiplies by `rfprm_scale`.
- `rfprm_scale` defaults to `(1.0u"W/m^2/K", 1e7u"Ω/m^2", 1e7u"Ω/m^2")`.
- `kwargs` is passed directly (as is) to the ODE `solve` call.

"""
function gen_sol_rf(KBB_log, po::ParamObjRF; rfprm_scale = rfprm_base_scale, kwargs...)
    newp = deepcopy(po)
    @reset newp.K_vwf = exp(KBB_log[1])*rfprm_scale[1]
    @reset newp.B_f = exp(KBB_log[2])*rfprm_scale[2]
    @reset newp.B_vw = exp(KBB_log[3])*rfprm_scale[3]
    newprob = ODEProblem(newp)
    sol = solve(newprob, Rodas3(); kwargs...)
    return sol
end
function gen_sol_rf(KBB_log, po::ParamObjRF, fitdat::PrimaryDryFit; rfprm_scale = rfprm_base_scale, kwargs...)
    newp = deepcopy(po)
    @reset newp.K_vwf = exp(KBB_log[1])*rfprm_scale[1]
    @reset newp.B_f = exp(KBB_log[2])*rfprm_scale[2]
    @reset newp.B_vw = exp(KBB_log[3])*rfprm_scale[3]
    newprob = ODEProblem(newp)
    sol = solve(newprob, Rodas3(); saveat=ustrip.(u"hr", fitdat.t), kwargs...)
    return sol
end

"""
    obj_KBB(KBB_log, pf; tweight=1.0, verbose=false)

Calculate the sum of squared error (objective function) for fitting parameters in the LC3 model.

# Arguments
- `KBB_log`: The nondimensional logs of `[Kvwf, Bf, Bvw]`; see [`LyoPronto.gen_sol_rf`](@ref).
- `pf`: Tuple of `(p::ParamObjRF, f::PrimaryDryFit)`
- `tweight`: (Optional) A weighting factor for the objective function. Default is `1.0`.
- `verbose`: (Optional) A boolean flag to enable verbose output. Default is `false`.
"""
function obj_KBB(KBB_log, pf; tweight=1.0, verbose=false)
    rtype = eltype(KBB_log)
    sol = gen_sol_rf(KBB_log, pf[1], pf[2])
    return obj_expT(sol, pf[2], tweight=tweight, verbose=verbose)::rtype
end

# @doc raw"""
#     obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u"hr", tweight=1, verbose=true)

# Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTTdat`.

# - `gen_sol` is a function taking [α, Kvwf, Bf, Bvw] and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
# - `tTTdat` is experimental temperature series, of the form `(time, Tf, Tvw)`, so with frozen and vial wall temperatures taken at the same time points. 
#     See also [`obj_tT_rf`](@ref) and [`obj_ttTT_rf`](@ref).
# - `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
# - `t_end` has a default value of `0.0u"hr"`, which (if left at default) is not included in the objective function.
# - `t_end` defaults to `missing`, in which case it is excluded from the objective.
#     If another value is passed, then the model drying time is compared to that value and included in the objective.
# """
# function obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=missing, tweight=1, verbose=true)
#     if any(fitprm .< 0)
#         return NaN
#     end
#     tloc, Tfdat, Tvwdat = tTTdat
#     pdfit = PrimaryDryFit(tloc, (Tfdat,), tloc, (Tvwdat,), t_end)
#     sol = gen_sol(fitprm)[1]
#     verbose && @info "gen_sol" fitprm 
#     return obj_expT(sol, pdfit; tweight=tweight, verbose=verbose)
# end


# @doc raw"""
#     obj_tT_rf(fitprm, gen_sol, tTdat; t_end=0.0u"hr", tweight=1, Tvw_end = 0.0u"K", verbose=true)

# Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTdat`.

# - `gen_sol` is a function taking `[α, Kvwf, Bf, Bvw]` and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
# - `tTTdat` is experimental temperature series, of the form `(time, Tf)`, so with frozen temperatures only. 
#     See also [`obj_tTT_rf`](@ref) and [`obj_ttTT_rf`](@ref).
# - `Tvw_end` defaults to `missing`, in which case vial wall temperatures are excluded from the objective. 
#     If another value is passed, then the final model vial wall temperature is compared to that value and included in the objective.
# - `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
# - `t_end` defaults to `missing`, in which case it is excluded from the objective.
#     If another value is passed, then the model drying time is compared to that value and included in the objective.
# """
# function obj_tT_rf(fitprm, gen_sol, tTdat; t_end=missing, tweight=1, Tvw_end = missing, verbose=true)
#     if any(fitprm .< 0)
#         return NaN
#     end
#     tloc, Tfdat = tTdat
#     pdfit = PrimaryDryFit(tloc, (Tfdat,), missing, Tvw_end, t_end)
#     sol = gen_sol(fitprm)[1]
#     verbose && @info "gen_sol" fitprm 
#     return obj_expT(sol, pdfit; tweight=tweight, verbose=verbose)
# end

# @doc raw"""
#     obj_ttTT_rf(fitprm, gen_sol, ttTTdat; t_end=missing, tweight=1, verbose=true)

# Evaluate an objective function which compares model solution with `fitprm` to experimental data in `tTdat`.

# - `gen_sol` is a function taking `[α, Kvwf, Bf, Bvw]` and returning a solution to lumped-capacitance microwave-assisted model; see [`gen_sol_rf_dim`](@ref).
# - `ttTTdat` is experimental temperature series, of the form `(time_Tf, time_vw, Tf, Tvw)`, so with Tf and Tvw having separate time points. 
#     This is useful if there is an early temperature rise in Tf, but Tvw continues to be reliable, so the model can fit to as much of Tvw as reasonable.
#     See also [`obj_tTT_rf`](@ref) and [`obj_tT_rf`](@ref).
# - `tweight` gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.
# - `t_end` defaults to `missing`, in which case it is excluded from the objective.
#     If another value is passed, then the model drying time is compared to that value and included in the objective.
# """
# function obj_ttTT_rf(fitprm, gen_sol, ttTTdat; t_end=missing, tweight=1, verbose=true)
#     if any(fitprm .< 0)
#         return NaN
#     end
#     tf, tvw, Tfdat, Tvwdat = ttTTdat
#     pdfit = PrimaryDryingFit(tf, (Tfdat,), tvw, (Tvwdat), t_end)
#     sol = gen_sol(fitprm)[1]
#     verbose && @info "gen_sol" fitprm 
#     return obj_expT(sol, pdfit; tweight=tweight, verbose=verbose)
# end
