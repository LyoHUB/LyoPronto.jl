
# It would be nice to use ConcreteStructs for this, but it struggles with my constructor.
struct PrimaryDryFit{T1, T2, T3, T4, T5, T6}
    t::T1
    Tfs::T2
    Tf_iend::T3# = [length(Tf) for Tf in Tfs]
    Tvws::T4# = missing
    Tvw_iend::T5# = (ismissing(Tvws) ? missing : [length(Tvw) for Tvw in Tvws])
    t_end::T6# = missing
    function PrimaryDryFit(t, Tfs, Tf_iend, Tvws, Tvw_iend, t_end)
        if !(t isa AbstractVector) || ~(first(t) isa Unitful.Time)
            throw(ArgumentError("t should be a vector of time points with units of time"))
        end
        if !(Tfs isa Tuple) || !(Tfs[1] isa AbstractVector) 
            throw(ArgumentError("Tfs should be a tuple of vectors of temperature measurements"))
        end
        if !ismissing(Tvws)
            if Tvws isa Number && !(Tvws isa Unitful.Temperature)
                throw(ArgumentError("If Tvws is a single value, it should be a temperature"))
            end
            if (Tvws isa AbstractVector) || ~(Tvws[1][1] isa Unitful.Temperature)
                throw(ArgumentError("If not a single value, Tvws should be a tuple of vectors of temperature measurements"))
            end
        end
        if !ismissing(t_end)
            if !(t_end isa Unitful.Time) && (t_end isa Tuple && !(first(t_end) isa Unitful.Time))
                throw(ArgumentError("t_end should be a time or a tuple of two times"))
            end
        end
        new{typeof.((t, Tfs, Tf_iend, Tvws, Tvw_iend, t_end))...}(
                t, Tfs, Tf_iend, Tvws, Tvw_iend, t_end)
    end
end
@doc """
PrimaryDryFit: a type for indicating how experimental data should be fit.

Provided constructors:

    PrimaryDryFit(t, Tfs, Tvws, t_end)
    PrimaryDryFit(t, Tfs; Tvws=missing, t_end=missing) = PrimaryDryFit(t, Tfs, Tvws, t_end)

Note that `Tvws` and `t_end` in the second constructors are keyword arguments, so either or both
can be left out.

The use of this struct is determined in large part by the implementation of 
[`LyoPronto.obj_expT`](@ref) and [`LyoPronto.err_expT!`](@ref). If a given field is not 
available, it will be set to `missing` and things should basically work. At least `t` and 
`Tfs` are expected to always be provided.

In the end, `Tfs` and `Tvws` are each stored as a tuple of vectors, but the constructors try to 
be flexible about allowing a single vector to be passed in place of a tuple of vectors.

The fields `Tf_iend` and `Tvw_iend` default to `[length(Tf) for Tf in Tfs]` and `[length(Tvw) for Tvw in Tvws]`, 
respectively, with one value for each temperature series; 
they are used to dictate if a given temperature series should
be truncated sooner than the full length in the fitting procedure.
This implies that all the temperature series correspond to the same
time points, then stop having measured values after a different number of measurements.
If a single value is given for `Tvws`, then it is taken to be an endpoint, and `Tvw_iend` will be `missing`.

`t_end` indicates an end of drying, particularly if taken from other measurements
(e.g. from Pirani-CM convergence). If set to `missing`, it is ignored in the
objective function. If set to a tuple of two times, then in the objective function any time 
in that window is not penalized; outside that window, squared error takes over, as for the 
single time case.

Principal Cases:
- Conventional: provide only `t, Tfs`
- Conventional with Pirani ending: provide `t, Tfs; t_end=...`
- RF with measured vial wall: provide `t, Tfs; Tvws=...`, 
- RF, matching model Tvw to experimental Tf[end] without measured vial wall: provide `t, Tfs, Tvws=...`
"""
PrimaryDryFit

# Primary constructor
function PrimaryDryFit(t, Tfs, Tvws, t_end) 
    if Tfs isa AbstractVector
        if eltype(Tfs) <: Number
            Tfs = (Tfs,)
        else
            Tfs = Tuple(Tfs...)
        end
    end
    if Tvws isa AbstractVector
        if eltype(Tvws) <: Number
            Tvws = (Tvws,)
        else
            Tvws = Tuple(Tvws...)
        end
    end
    if t_end isa Tuple
        t_end = extrema(t_end)
    end
    PrimaryDryFit(t, Tfs, [length(Tf) for Tf in Tfs], Tvws, 
    ((ismissing(Tvws) || Tvws isa Number) ? missing : [length(Tvw) for Tvw in Tvws]),
    t_end)
end
# Convenience constructors
PrimaryDryFit(t, Tfs; Tvws=missing, t_end=missing) = PrimaryDryFit(t, Tfs, Tvws, t_end)

function Base.:(==)(p1::PrimaryDryFit, p2::PrimaryDryFit)
    cond1 = p1.t == p2.t
    cond2 = p1.Tfs == p2.Tfs
    cond3 = p1.Tf_iend == p2.Tf_iend
    cond4 = ismissing(p1.Tvws) ? ismissing(p2.Tvws) : (p1.Tvws == p2.Tvws)
    cond5 = ismissing(p1.Tvw_iend) ? ismissing(p2.Tvw_iend) : (p1.Tvw_iend == p2.Tvw_iend)
    cond6 = ismissing(p1.t_end) ? ismissing(p2.t_end) : (p1.t_end == p2.t_end)
    return all([cond1, cond2, cond3, cond4, cond5, cond6])
end


"""
    $(SIGNATURES)

Evaluate an objective function which compares model solution computed by `sol` to experimental data in `pdfit`.

- `sol` is a solution to an appropriate model; see [`gen_sol_pd`](@ref) for a helper.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
- `Tvw_weight = 1` gives the weighting of Tvw in the objective, as compared to Tf.

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for the lumped capacitance model (see [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `pdfit`, squared error is computed for each separately then summed; likewise for `Tvw`.

I've considered writing several methods and dispatching on `pdfit` somehow, which would be cool and might individually be easier to read. But control flow might be harder to document and explain, and this should work just fine.
"""
function obj_expT(sol::ODESolution, pdfit::PrimaryDryFit; 
    tweight=1.0, verbose = false, Tvw_weight=1.0)
    if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
        verbose && @warn "ODE solve did not reach end of drying. Either parameters are bad, or tspan is not large enough." sol.retcode sol.prob.p.hf0 sol[end]
        return Inf
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

obj_expT(sol::Val{NaN}, pdfit; kwargs...) = Inf
function obj_expT(sol, pdfit; verbose=false, kwargs...) 
    verbose && @warn "`obj_expT` got passed improper args. Might not be a problem, but check." sol
    # In some cases, inputs are so bad it's not worth an ODE solve, so this method
    # provides an escape hatch for NaN returns instead of crashing.
    if !isfinite(sol) 
        return Inf
    end
    error("Improper call to `obj_expT`.")
end


"""
    $(SIGNATURES)

Compute the number of data points available in `pdfit` for comparison to model solution.

This is useful for caching a residual vector for least-squares fitting, e.g. with [`err_expT!`](@ref) and [`nls_pd!`](@ref).
"""
function num_errs(pdfit)
    # Count the number of errors in the PrimaryDryFit object
    Tvw_len = ismissing(pdfit.Tvws) ? 0 : (ismissing(pdfit.Tvw_iend) ? 1 : sum(pdfit.Tvw_iend))
    nerr = sum(pdfit.Tf_iend) + Tvw_len + (ismissing(pdfit.t_end) ? 0 : 1)
    return nerr
end

const errexpT_doc = """
Evaluate the error between model solution `sol` and experimental data in `pdfit`.

The in-place version `err_expT!` fills the `errs` vector with the errors, while the non-in-place version `err_expT` returns a new vector of errors.
In-place `err_expT!` thus doesn't allocate, but requires `errs` to have length `num_errs(pdfit)`.

In contrast to `obj_expT()`, which sums all the squared residuals, this function fills the passed array
`errs` with each separate residual, which is suited for least squares algorithms.
- `errs` is a vector of length `num_errs(pdfit)`, which this function fills with the errors.
-`sol` is a solution to an appropriate model; see [`gen_sol_pd`](@ref) for a helper function.
- `pdfit` is an instance of `PrimaryDryFit`, which contains some information about what to compare.
- `tweight = 1` gives the weighting (in K^2/hr^2) of the total drying time in the objective, as compared to the temperature error.
Each time series, plus the end time, is given equal weight by dividing by its length; error is given in K (but `ustrip`ped).

Note that if `pdfit` has vial wall temperatures (i.e. `ismissing(pdfit.Tvws) == false`), the third-index variable in `sol` is assumed to be temperature, as is true for solutions with [`ParamObjRF`](@ref).

If there are multiple series of `Tf` in `pdfit`, error is computed for each separately; likewise for `Tvw`.
"""

"""
    $(SIGNATURES)

$errexpT_doc
"""
function err_expT!(errs, sol::ODESolution, pdfit; tweight=1, verbose = false)
    if length(errs) != num_errs(pdfit)
        error("Wrong length of cached residual vector.")
    end
    # errs .= 0.0 # If indexing is handled correctly, this should not be necessary.
    not_avail_err = 0.0 # an error value to return for points where the solution is unavailable, e.g. if the model dries faster
    if sol.retcode != ReturnCode.Terminated || length(sol.u) <= 2
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
        errs .= Inf
        return
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
        errs[last_ind+1:last_ind+i_solstart] .= not_avail_err
        errs[last_ind+i_solstart:last_ind+trim] .= ustrip.(u"K", Tferrs)
        errs[last_ind+trim+1:last_ind+itf] .= not_avail_err
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
                errs[last_ind+1:last_ind+i_solstart] .= not_avail_err
                errs[last_ind+i_solstart:last_ind+trim] .= ustrip.(u"K", Tvw_errs)
                errs[last_ind+trim+1:last_ind+itvw] .= not_avail_err
                last_ind += itvw
            end
        end
    end

    # Add end time at end of array, if present
    if !ismissing(pdfit.t_end)
        if pdfit.t_end isa Tuple # See if is inside window and scale appropriately
            mid_t = (pdfit.t_end[1] + pdfit.t_end[2]) / 2.0
            if tmd < pdfit.t_end[1]
                t_err = (mid_t - tmd)
            elseif tmd > pdfit.t_end[2]
                t_err = (mid_t - tmd)
            else # Inside window, so no error
                t_err = 0.0u"hr"
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
err_expT!(errs, sol::Float64, pdfit; kwargs...) = isnan(sol) ? errs .= Inf : error("Unexpected state")
err_expT!(errs, sol::Val{NaN}, pdfit; kwargs...) = errs .= Inf;

"""
    $(SIGNATURES)

$errexpT_doc
"""
function err_expT(sol, pdfit; tweight=1, verbose = false)
    if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
        verbose && @info "ODE solve failed or incomplete, probably." sol.retcode sol[1, :]
        return [Inf]
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
    # Compute temperature errors for all frozen temperatures
    # for (Tf, itf) in zip(pdfit.Tfs, pdfit.Tf_iend)
    errs = mapreduce(vcat, pdfit.Tfs, pdfit.Tf_iend) do Tf, itf
        trim = min(itf, length(Tfmd))
        Tferrs = (Tf[i_solstart:trim] .- Tfmd[begin:trim-i_solstart+1])/sqrt(trim-i_solstart+1)
        # append!(errs, ustrip.(u"K", Tferrs))
        return ustrip.(u"K", Tferrs)
    end

    # If present, vcat vial wall temperatures
    if !ismissing(pdfit.Tvws) # At least one vial wall temperature
        if ismissing(pdfit.Tvw_iend) # Only an endpoint temperature provided
            Tvwend = pdfit.Tvws
            Tvw_err = sol[3, end]*u"K" - uconvert(u"K", Tvwend)
            push!(errs, ustrip(u"K", Tvw_err))
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
                append!(errs, ustrip.(u"K", Tvw_errs))
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
                t_err = 0.0u"hr"
            end
        else
            t_err = (pdfit.t_end - tmd)
        end
        push!(errs, ustrip(u"hr", t_err*tweight))
    end
    verbose && @info "loss call" tmd size(errs) sum(abs2.(errs))
    return errs
end
