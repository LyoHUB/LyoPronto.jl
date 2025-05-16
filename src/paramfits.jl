export gen_sol_pd, obj_pd, gen_nsol_pd, objn_pd
export KRp_transform_basic, K_transform_basic, Rp_transform_basic
export KBB_transform_basic
export obj_expT 
export err_expT, num_errs, nls_pd

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
sol = solve(prob, Rodas4(autodiff=AutoForwardDiff(chunksize=2)); saveat, kwargs...)
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
    sol = solve(prob, odealg_chunk2; saveat, kwargs...)
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

Calculate the sum of squared error (objective function) for fitting parameters to primary drying data.
This directly calls [`gen_sol_pd`](@ref), then [`obj_expT`](@ref), so see those docstrings.
"""
function nls_pd(errs, fitlog, tpf; tweight=1.0, verbose=false)
    sol = gen_sol_pd(fitlog, tpf...)
    return err_expT(errs, sol, tpf[3]; tweight, verbose)
end

"""
    $(SIGNATURES)

Evaluate an objective function which compares model solution computed by `sol` to experimental data in `pdfit`.

- `sol` is a solution to an appropriate model; see [`gen_sol_pd`](@ref) for a helper.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
- `Tvw_weight = 1` gives the weighting of Tvw in the objective, as compared to Tf.

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for the lumped capacitance model (see [ParamObjRF`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.

I've considered writing several methods and dispatching on `pdfit` somehow, which would be cool and might individually be easier to read. But control flow might be harder to document and explain, and this should work just fine.
"""
function obj_expT(sol::ODESolution, pdfit::PrimaryDryFit; 
    tweight=1.0, verbose = false, Tvw_weight=1.0)
    if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
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
    for (Tf, iend) in zip(pdfit.Tfs, pdfit.Tf_iend)
        trim = min(iend, length(Tfmd))
        Tfobj += sum(abs2, (Tf[i_solstart:trim] .- Tfmd[begin:trim-i_solstart+1]))/(trim-i_solstart+1)
    end
    if ismissing(pdfit.Tvws) # No vial wall temperatures
        Tvwobj = 0.0u"K^2"
    elseif ismissing(pdfit.Tvw_iend) # Only an endpoint temperature provided
        Tvwend = pdfit.Tvws
        Tvwobj = (sol[3, end]*u"K" - uconvert(u"K", Tvwend))^2
    else # Regular case of fitting to at least one full temperature series
        if preinterp
            Tvwmd = sol[3, begin:end-1].*u"K" # Leave off last time point because is end time
        else
            vwtrim = sol.t[begin]*u"hr" .< pdfit.t .< tmd
            tvw_trim = pdfit.t[vwtrim]
            Tvwmd = sol.(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
        end
        # Compute temperature objective for all vial wall temperatures
        Tvwobj = 0.0u"K^2"
        for (Tvw, iend) in zip(pdfit.Tvws, pdfit.Tvw_iend)
            trim = min(iend, length(Tvwmd))
            Tvwobj += sum(abs2, (Tvw[i_solstart:trim] .- Tvwmd[begin:trim-i_solstart+1]))/(trim-i_solstart+1)
        end
    end
    if ismissing(pdfit.t_end) # No drying time provided
        tobj = 0.0u"hr^2"
    elseif pdfit.t_end isa Tuple # See if is inside window and scale appropriately
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
    # In some cases, inputs are so bad it's not worth an ODE solve, so this method
    # provides an escape hatch for NaN returns instead of crashing.
    if isnan(sol) 
        return NaN
    end
    error("Improper call to `obj_expT`.")
end

function num_errs(pdfit)
    # Count the number of errors in the PrimaryDryFit object
    Tvw_len = ismissing(pdfit.Tvws) ? 0 : (ismissing(pdfit.Tvw_iend) ? 1 : sum(pdfit.Tvw_iend))
    nerr = sum(pdfit.Tf_iend) + Tvw_len + (ismissing(pdfit.t_end) ? 0 : 1)
    return nerr
end
"""
    $(SIGNATURES)

Evaluate the error between model solution `sol` to experimental data in `pdfit`.

In contrast to `obj_expT()`, this function makes an array of all the errors, which would be squared and summed to produce an objective function.
- `errs` is a vector of length `num_errs(pdfit)`, which this function fills with the errors.
-`sol` is a solution to an appropriate model; see [`gen_sol_pd`](@ref) for a helper function.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
Each time series, plus the end time, is given equal weight by dividing by its length; error is given in K (but `ustrip`ped).

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for solutions with [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `pdfit`, error is computed for each separately; likewise for `Tvw`.
"""
function err_expT(errs, sol, pdfit; tweight=1, verbose = false)
    if length(errs) != num_errs(pdfit)
        error("Wrong length of provided error vector.")
    end
    # errs .= 0.0 # If indexing is handled correctly, this should not be necessary.
    if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
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
    # Initialize error array with frozen temperatures
    # Compute temperature errors for all frozen temperatures
    # errs = Float64[]
    last_ind = 0
    for (Tf, itf) in zip(pdfit.Tfs, pdfit.Tf_iend)
    # errs = mapreduce(vcat, pdfit.Tfs, pdfit.Tf_iend) do Tf, itf
        trim = min(itf, length(Tfmd))
        Tferrs = (Tf[i_solstart:trim] .- Tfmd[begin:trim-i_solstart+1])/sqrt(trim-i_solstart+1)
        errs[last_ind+1:last_ind+i_solstart] .= 0.0
        errs[last_ind+i_solstart:last_ind+trim] .= ustrip.(u"K", Tferrs)
        errs[last_ind+trim+1:last_ind+itf] .= 0.0
        last_ind += itf
    end

    # If present, vcat vial wall temperatures
    if !ismissing(pdfit.Tvws) # At least one vial wall temperature
        if ismissing(pdfit.Tvw_iend) # Only an endpoint temperature provided
            Tvwend = pdfit.Tvws
            Tvw_err = sol[3, end]*u"K" - uconvert(u"K", Tvwend)
            errs[last_ind+1] = ustrip(u"K", Tvw_err)
            last_ind += 1
        else # Regular case of fitting to at least one full temperature series
            if preinterp
                Tvwmd = sol[3, begin:end-1].*u"K" # Leave off last time point because is end time
            else
                vwtrim = sol.t[begin]*u"hr" .< pdfit.t .< tmd
                tvw_trim = pdfit.t[vwtrim]
                Tvwmd = sol.(ustrip.(u"hr", tvw_trim), idxs=3).*u"K"# .- 273.15
            end
            # Compute temperature objective for all vial wall temperatures
            for (Tvw, itvw) in zip(pdfit.Tvws, pdfit.Tvw_iend) 
                trim = min(itvw, length(Tvwmd))
                Tvw_errs = (Tvw[i_solstart:trim] .- Tvwmd[begin:trim-i_solstart+1])/sqrt(trim-i_solstart+1)
                errs[last_ind+1:last_ind+i_solstart] .= 0.0
                errs[last_ind+i_solstart:last_ind+trim] .= ustrip.(u"K", Tvw_errs)
                errs[last_ind+trim+1:last_ind+itvw] .= 0.0
                last_ind += itvw
            end
        end
    end

    # Concatenate end time to array, if present
    if !ismissing(pdfit.t_end)
        if pdfit.t_end isa Tuple # See if is inside window and scale appropriately
            mid_t = (pdfit.t_end[1] + pdfit.t_end[2]) / 2.0
            if tmd < pdfit.t_end[1]
                t_err = (mid_t - tmd)
            elseif tmd > pdfit.t_end[2]
                t_err = (mid_t - tmd)
            else # Inside window, so no error
                t_err = 0.0u"hr^2"
            end
        else
            t_err = (pdfit.t_end - tmd)
        end
        errs[last_ind+1] = ustrip(u"hr", t_err*tweight)
        if last_ind + 1 != length(errs)
            error("Indexing problems...")
        end
    end
    verbose && @info "loss call" tmd size(errs) sum(abs2.(errs))
    return nothing
end
function err_expT(sol::ODESolution, pdfit::PrimaryDryFit; kwargs...)
    errs = fill(0.0, num_errs(pdfit))
    err_expT(errs, sol, pdfit; kwargs...)
    return errs
end
