

"""
    $(SIGNATURES)

Construct a typical transform for fitting both Kshf and Rp.
"""
function KRp_transform_basic(Kshfg, R0g, A1g, A2g)
    t1 = K_transform_basic(Kshfg)
    t2 = Rp_transform_basic(R0g, A1g, A2g)
    return merge(t1, t2)
end
"""
    $(SIGNATURES)

Construct a typical transform for fitting Rp.
"""
function Rp_transform_basic(R0g, A1g, A2g)
    tr = as((;
        Rp = as(RpFormFit, as((;
            R0 = TVScale(R0g) ∘ TVExp(),
            A1 = TVScale(A1g) ∘ TVExp(),
            A2 = TVScale(A2g) ∘ TVExp(),
            )),)
        ))
    return tr
end
"""
    $(SIGNATURES)

Construct a typical transform for fitting Kshf (a.k.a. Kv).
"""
function K_transform_basic(Kshfg)
    tr = as((;Kshf = as(ConstPhysProp, (TVScale(Kshfg) ∘ TVExp(),))))
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
Construct a bounded transform for fitting Kvwf, Bf, and Bvw (as for a microwave cycle).

`Kvwf_scalefac`, `Bf_scalefac`, and `Bvw_scalefac` are used to provide upper and lower 
bounds on the fitted parameter, as `(Kvwfg/Kvwf_scalefac, Kvwfg*Kvwf_scalefac)`, etc.
This is enforced with a logistic transform scaled and shifted appropriately.
"""
function KBB_transform_bounded(Kvwfg, Bfg, Bvwg; Kvwf_scalefac=1e2, Bf_scalefac=1e4, Bvw_scalefac=1e4)
    tr = as((Kvwf = TVScale(Kvwfg*Kvwf_scalefac) ∘ TVLogistic() ∘ TVShift(logit(inv(Kvwf_scalefac))),
        Bf = TVScale(Bfg*Bf_scalefac) ∘ TVLogistic() ∘ TVShift(logit(inv(Bf_scalefac))),
        Bvw = TVScale(Bvwg*Bvw_scalefac) ∘ TVLogistic() ∘ TVShift(logit(inv(Bvw_scalefac))) ))
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
!isnothing(badprms) && badprms(new_params) && return Val(NaN)
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
function gen_sol_pd(fitlog, tr, po; saveat=[], badprms=nothing, kwargs...)
    fitprm = transform(tr, fitlog)
    prms = setproperties(po, fitprm)
    !isnothing(badprms) && badprms(prms) && return Val(NaN)
    prob = ODEProblem(prms; tspan=(0.0, 1000.0))
    sol = solve(prob, odealg_chunk2; saveat, kwargs...)
    return sol
end
"$(SIGNATURES)"
function gen_sol_pd(fitlog, tr, po, fitdat; badprms=nothing, kwargs...)
    sol = gen_sol_pd(fitlog, tr, po; saveat=ustrip.(u"hr", fitdat.t), badprms, kwargs...)
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
function gen_nsol_pd(fitlog, tr, pos, fitdats; badprms=nothing, kwargs...)
    saveats = [ustrip.(u"hr", fitdat.t) for fitdat in fitdats]
    return gen_nsol_pd(fitlog, tr, pos; saveats, badprms, kwargs...)
end
function gen_nsol_pd(fitlog, tr, pos; saveats=fill([], length(pos)), badprms=nothing, kwargs...)
    fitprm = transform(tr, fitlog)
    if pos isa ParamObj # only one param object...
        pos = repeat([pos], length(saveats))
    end
    if hasproperty(fitprm, :separate) && hasproperty(fitprm, :shared)
        sep_inds = get(fitprm, :sep_inds, 1:length(fitprm.separate))
        if length(sep_inds) != length(saveats)
            error("Length of either the transformed variable or `sep_inds` does not match fitdats.")
        end
        prms =  [setproperties(po, merge(s, fitprm.shared)) for (po, s) in zip(pos, fitprm.separate[sep_inds])]
    else
        prms = [setproperties(po, fitprm) for po in pos]
    end
    if length(prms) != length(pos)
        error("Length of transformed variable and pos do not match.")
    end
    if !isnothing(badprms) && any([badprms(p) for p in prms])
       return fill(Val(NaN), length(prms))
    end
    sols = map(prms, saveats) do prm, saveat
        prob = ODEProblem(prm; tspan=(0.0, 1000.0))
        return solve(prob, odealg_chunk2; saveat, kwargs...)
    end
    return sols
end

"""
    $(SIGNATURES)

Calculate the sum of squared error (objective function) for fitting parameters to primary drying data.
This directly calls [`gen_sol_pd`](@ref), then [`obj_expT`](@ref), so see those docstrings.
"""
function obj_pd(fitlog, tpf; tweight=1.0, Tvw_weight=1.0, badprms=nothing, verbose=false)
    sol = gen_sol_pd(fitlog, tpf...; badprms)
    return obj_expT(sol, tpf[3]; tweight, Tvw_weight, verbose)
end

"""
    $(SIGNATURES)

Calculate the sum of squared error (objective function) for fitting parameters to primary drying data.
This directly calls [`gen_nsol_pd`](@ref), then [`obj_expT`](@ref), so see those docstrings.
"""
function objn_pd(fitlog, tpf; tweight=1.0, Tvw_weight=1.0, badprms=nothing, verbose=false)
    sols = gen_nsol_pd(fitlog, tpf...; badprms)
    obj = mapreduce(+, sols, tpf[3]) do sol, fitdat
        obj_expT(sol, fitdat; tweight, Tvw_weight, verbose)
    end
    return obj
end

"""
    $(SIGNATURES)

Calculate the errors for fitting parameters to primary drying data.
This directly calls [`gen_sol_pd`](@ref), then [`err_expT!`](@ref), so see those docstrings.
"""
function nls_pd!(errs, fitlog, tpf; tweight=1.0, verbose=false)
    sol = gen_sol_pd(fitlog, tpf...)
    return err_expT!(errs, sol, tpf[3]; tweight, verbose)
end
"""
    $(SIGNATURES)

Calculate the errors for fitting parameters to primary drying data.
This directly calls [`gen_sol_pd`](@ref), then [`err_expT`](@ref), so see those docstrings.
"""
function nls_pd(fitlog, tpf; tweight=1.0, verbose=false)
    sol = gen_sol_pd(fitlog, tpf...)
    return err_expT(sol, tpf[3]; tweight, verbose)
end

# Prepare a fitting nonlinear function with some sensible defaults
function NonlinearFunction(fitdat::PrimaryDryFit; tweight=1.0, verbose=false)
    if tweight == 1.0 && verbose == false
        NonlinearFunction{true, SciMLBase.FullSpecialize}(nls_pd!, 
            resid_prototype=zeros(num_errs(fitdat)))
    else
        f = (e, f, tpf) -> nls_pd!(e, f, tpf; tweight, verbose)
        NonlinearFunction{true, SciMLBase.FullSpecialize}(f,
            resid_prototype=zeros(num_errs(fitdat)))
    end
end
