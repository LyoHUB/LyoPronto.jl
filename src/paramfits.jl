export gen_sol_pd, obj_pd, gen_nsol_pd, objn_pd
export KRp_transform_basic, K_transform_basic, Rp_transform_basic
export KBB_transform_basic
export obj_expT, err_expT 

"""
    $(SIGNATURES)

Construct a typical transform for fitting both Kshf and Rp.
"""
function KRp_transform_basic(Kshfg, R0g, A1g, A2g)
    t1 = K_transform_basic(Kshfg)
    t2 = Rp_transform_basic(R0g, A1g, A2g)
    return as(merge(t1.transformations, t2.transformations))
end
"""
    $(SIGNATURES)

Construct a typical transform for fitting Rp.
"""
function Rp_transform_basic(R0g, A1g, A2g)
    tr = as((
        Rp = as((
            R0 = TVScale(R0g) ∘ TVExp(),
            A1 = TVScale(A1g) ∘ TVExp(),
            A2 = TVScale(A2g) ∘ TVExp(),
            )),
        ))
    return tr
end
"""
    $(SIGNATURES)

Construct a typical transform for fitting Kshf (a.k.a. Kv).
"""
function K_transform_basic(Kshfg)
    tr = as((Kshf = ConstWrapTV() ∘ TVScale(Kshfg) ∘ TVExp(),))
    return tr
end
"""
    $(SIGNATURES)
Construct a typical transform for fitting Kvwf, Bf, and Bvw (as for a microwave cycle).
"""
function KBB_transform_basic(Kvwfg, Bfg, Bvwg)
    tr = as((Kvwf = TVScale(Kvwfg) ∘ TVExp(), 
        Bf = TVScale(Bfg) ∘ TVExp(),
        Bvw = TVScale(Bvwg) ∘ TVExp(),))
    return tr
end



"""
    $(SIGNATURES)

Simulate primary drying, given a vector of parameter guesses, a mapping `tr` from `fitlog` to named coefficients, and other parameters in `po`.

The equations used are determined by the type of `po`, which (with the magic of dispatch)
is used to set up an ODE system.

`tr` should be a `TransformTuple` object, from TransformVariables, which maps e.g. a vector of 3
real numbers to a NamedTuple with `R0, A1, A2` as keys and appropriate Unitful dimensions on the values.
This small function runs 
```
fitprm = transform(tr, fitlog)
new_params = setproperties(po, fitprm)
prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
sol = solve(prob, Rodas3(); saveat, kwargs...)
```
which is wrapped to avoid code duplication.

So, to choose which parameters to fitting, all that is necessary is to provide an appropriate transform `tr`
and add a method of `setproperties` for the desired parameters.
Therefore this function can be used for both K-Rp fitting, or just Rp, or just a subset of the 3 Rp coefficients.

If given, `fitdat` is used to set `saveat` for the ODE solution.

Other `kwargs` are passed directly (as is) to the ODE `solve` call.
"""
function gen_sol_pd(fitlog, tr, po; saveat=[], kwargs...)
    fitprm = transform(tr, fitlog)
    prms = setproperties(po, fitprm)
    prob = ODEProblem(prms; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3(); saveat, kwargs...)
    return sol
end
"$(SIGNATURES)"
function gen_sol_pd(fitlog, tr, po, fitdat; kwargs...)
    sol = gen_sol_pd(fitlog, tr, po; saveat=ustrip.(u"hr", fitdat.t), kwargs...)
    return sol
end

"""
    $(SIGNATURES)

Generate multiple solutions at once.

If the transformation `tr` makes something (e.g. NamedTuple) with properties `separate` and 
`shared`, then one each of `separate` is combined with `shared`, then they are matched up with each
element of `pos` and `fitdats`. If `pos` is a single object, it is repeated.
Further, if `tr` also has a field `sep_inds`, then those indices are used to map `separate` 
to the sets of `pos` and `fitdats`. This is useful for when e.g. 2 sets of `separate` parameters 
are to be applied across 5 different experiments.

"""
function gen_nsol_pd(fitlog, tr, pos, fitdats; kwargs...)
    saveats = [ustrip.(u"hr", fitdat.t) for fitdat in fitdats]
    return gen_nsol_pd(fitlog, tr, pos; saveats, kwargs...)
end
function gen_nsol_pd(fitlog, tr, pos; saveats=fill([], length(pos)), kwargs...)
    fitprm = transform(tr, fitlog)
    if pos isa ParamObj # only one param object...
        pos = repeat([pos], length(fitdats))
    end
    if hasproperty(fitprm, :separate) && hasproperty(fitprm, :shared)
        # if length(fitprm.separate) == 0
        #     prms =  [setproperties(po, fitprm.shared) for po in pos, fitprm.separate)]
        sep_inds = hasproperty(fitprm, :sep_inds) ? fitprm.sep_inds : 1:length(fitprm.separate)
        if length(sep_inds) != length(saveats)
            error("Length of either the transformed variable or `sep_inds` does not match fitdats.")
        end
            # prms =  [setproperties(po, merge(s, fitprm.shared)) for (po, s) in zip(pos, fitprm.separate)]
            # prms =  [setproperties(po, merge(s, fitprm.shared)) for (po, s) in zip(pos, fitprm.separate)]
        prms =  [setproperties(po, merge(s, fitprm.shared)) for (po, s) in zip(pos, fitprm.separate[sep_inds])]
    else
        prms = [setproperties(po, fitprm) for po in pos]
    end
    if length(prms) != length(pos)
        error("Length of transformed variable and pos do not match.")
    end
    # Allow for one shared parameter fit to several experiments at once
    # prms = setproperties(pos, fitprm)
    sols = map(prms, saveats) do prm, saveat
        prob = ODEProblem(prm; tspan=(0.0, 1000.0))
        soli = solve(prob, odealg_chunk2; saveat, kwargs...)
    end
    return sols
end

# function setproperties(pos::Vector{ParamObj}, patch::NamedTuple)
#     if hasproperty(patch, :separate) && hasproperty(patch, :shared)
#         return [setproperties(po, merge(s, patch.shared)) for (po, s) in zip(pos, patch.separate)]
#     else
#         return [setproperties(po, patch) for po in pos]
#     end
# end
# function setproperties(pos, patch::NamedTuple{(:separate, :shared)})
#     if length(pos) != length(patch.separate)
#         error("Length of pos and patch.separate do not match.")
#     end
#     map(pos, patch.separate) do po, s
#         setproperties(po, merge(s, patch.shared))
#     end
#     # return [setproperties(po, merge(s, patch.shared)) for (po, s) in zip(pos, patch.separate)]
# end

"""
    $(SIGNATURES)

Calculate the sum of squared error (objective function) for fitting parameters to primary drying data.
This directly calls [`gen_sol_pd`](@ref), then [`obj_expT`](@ref), so see those docstrings.
"""
function obj_pd(fitlog, tpf; tweight=1.0, verbose=false)
    sol = gen_sol_pd(fitlog, tpf...)
    return obj_expT(sol, tpf[3]; tweight, verbose)
end

function objn_pd(fitlog, tpf; tweight=1.0, verbose=false)
    sols = gen_nsol_pd(fitlog, tpf...)
    obj = mapreduce(+, sols, tpf[3]) do sol, fitdat
        obj_expT(sol, fitdat; tweight, verbose)
    end
    return obj
end

"""
    $(SIGNATURES)

Evaluate an objective function which compares model solution computed by `sol` to experimental data in `pdfit`.

- `sol` is a solution to an appropriate model; see [`gen_sol_Rp`](@ref), [`gen_sol_KRp`](@ref), and [`gen_sol for some helper functions for this.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
- `Tvw_weight = 1` gives the weighting of Tvw in the objective, as compared to Tf.

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for [`gen_sol_rf`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.

I've considered writing several methods and dispatching on `pdfit` somehow, which would be cool and might individually be easier to read. But control flow might be harder to document and explain, and this should work just fine.
"""
function obj_expT(sol::ODESolution, pdfit::PrimaryDryFit{TT1, TT2, TT3, TT4, TT5, TTvw, TTvwi, Tte};
     tweight=1.0, verbose = false, Tvw_weight=1.0) where {TT1, TT2, TT3, TT4, TT5, TTvw, TTvwi, Tte}
    if sol.retcode !== ReturnCode.Terminated
        # hf_end = sol[1, end]*u"cm"
        verbose && @warn "ODE solve did not reach end of drying. Either parameters are bad, or tspan is not large enough." sol.retcode sol.prob.p.hf0 sol[end]
        return NaN
    end
    tmd = sol.t[end].*u"hr"
    nt = length(sol.t) - 1
    i_solstart = searchsortedfirst(pdfit.t, sol.t[begin]*u"hr") 
    # Identify if the solution is pre-interpolated to the time points in pdfit.t
    preinterp = true
    for i in 1:nt
        if ~(sol.t[i]*u"hr" ≈ pdfit.t[i_solstart + i - 1])
            preinterp = false
            break
        end
    end
    # @info "here" preinterp i_solstart sol.t[begin:begin+4] pdfit.t[i_solstart:i_solstart+4]

    # Compute temperature objective for all frozen temperatures
    if preinterp
        Tfmd = sol[2, begin:end-1].*u"K" # Leave off last time point because is end time
    else
        ftrim = sol.t[begin]*u"hr" .< pdfit.t .< tmd
        tf_trim = pdfit.t[ftrim]
        Tfmd = sol.(ustrip.(u"hr", tf_trim), idxs=2).*u"K"
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
        Tfobj += sum(abs2, (pdfit.Tfs[j][i_solstart:trim] - Tfmd[begin:trim-i_solstart+1]))/(trim-i_solstart+1)
    end
    if TTvw == Missing # No vial wall temperatures, encoded in type
        Tvwobj = 0.0u"K^2"
    elseif TTvwi == Missing # Endpoint only temperature, encoded in type
        Tvwend = pdfit.Tvws
        Tvwobj = (sol[3, end]*u"K" - uconvert(u"K", Tvwend))^2
    else # Regular case of fitting to at least one full temperature series
        if preinterp
            Tvwmd = sol[3, begin:end-1].*u"K" # Leave off last time point because is end time
        else
            vwtrim = sol.t[begin]*u"hr" .< pdfit.t .< tmd
            tvw_trim = pdfit.t[vwtrim]
            Tvwmd = sol(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
        end
        # Compute temperature objective for all vial wall temperatures
        Tvwobj = 0.0u"K^2"
        for (j, iend) in enumerate(pdfit.Tvw_iend)
            trim = min(iend, length(Tvwmd))
            Tvwobj += sum(abs2, (pdfit.Tvws[j][i_solstart:trim] .- Tvwmd[begin:trim-i_solstart+1]))/(trim-i_solstart+1)
        end
    end
    if Tte == Missing # No drying time provided: encoded in type
        tobj = 0.0u"hr^2"
    elseif Tte <: Tuple # See if is inside window and scale appropriately
        mid_t = (pdfit.t_end[1] + pdfit.t_end[2]) / 2.0
        if tmd < pdfit.t_end[1]
            tobj = (mid_t - tmd)^2
        elseif tmd > pdfit.t_end[2]
            tobj = (mid_t - tmd)^2
        else # Inside window, so no error
            tobj = 0.0u"hr^2"
        end
    else # Compare to a single drying time
        tobj = (pdfit.t_end - tmd)^2
    end
    verbose && @info "loss call" tmd tobj Tfobj Tvwobj 
    return ustrip(u"K^2", Tfobj + Tvw_weight*Tvwobj) + tweight*ustrip(u"hr^2", tobj)
end

function obj_expT(sol, pdfit; verbose=false, kwargs...)
    verbose && @warn "`obj_expT` got passed improper args. Might not be a problem, but check." sol
    if isnan(sol)
        return NaN
    end
    error("Improper call to `obj_expT`.")
end

"""
    $(SIGNATURES)

Evaluate the error between model solution `sol` to experimental data in `pdfit`.

In contrast to `obj_expT()`, this function returns an array of all the errors, which would be squared and summed to produce an objective function.

- `sol` is a solution to an appropriate model; see [`gen_sol_Rp`](@ref), [`gen_sol_KRp`](@ref), and [`gen_sol_rf`](@ref) for some helper functions for this.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
Each time series, plus the end time, is given equal weight by dividing by its length; error is given in K (but `ustrip`ped).

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for [`gen_sol_rf`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.
"""
function err_expT(sol, pdfit; verbose = false)
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