
function new_Tsh_ramp(orig, new_setpt)
    return RampedVariable([orig.setpts[1], new_setpt|>u"K"], orig.ramprates[1])
end

"""
    $(SIGNATURES)

Calculate a single shelf temperature isotherm for a range of pressures.

Returns a table of results, with each row corresponding to a different pch in pch_range. Each row contains the following fields:
- `po`: The parameter object for the specific case
- `max_md`: The maximum mass flow from the solution
- `ave_md`: The average mass flow from the solution
- `end_md`: The mass flow at the end of primary drying
- `max_Tf`: The maximum product temperature from the solution
- `dry_time`: The total time to reach the end of primary drying
"""
function ds_shelf_isotherm(Tsh_iso, pch_range, po)
    tab = Table(map(pch_range) do pch_i
        # Set up the specific case and simulate
        Tsh = new_Tsh_ramp(po.Tsh, Tsh_iso)
        pch = RampedVariable(pch_i)
        po_i = setproperties(po, (;Tsh, pch))
        sol = solve(ODEProblem(po_i), odealg_chunk2)
        # Compute relevant quantities from the solution
        summary = summary_md_Q(sol)
        max_md = maximum(summary.md)
        dry_time = summary.t[end]
        # The total mass loss can likely be computed directly from Vfill and rho_solution,
        # but this is more general in case porosity gets defined differently, etc.
        del_t = diff(summary.t)
        weight_t = (vcat(0.0u"hr", del_t) + vcat(del_t, 0.0u"hr"))/2
        total_m = sum(summary.md .* weight_t)
        m = po.hf0*po.Ap * (po_i.ρsolution - po_i.csolid)
        ave_md = total_m/dry_time
        end_md = summary.md[end]
        max_Tf = maximum(sol[2,:])*u"K"
        return (;po=po_i, pch=pch_i, Tsh=Tsh_iso, max_md, ave_md, end_md, max_Tf, dry_time)
    end)
    return tab
end

# First, a generic fallback which doesn't require model manipulation
function ds_product_isotherm(Tpr, pch_range, po)
    tab = Table(map(pch_range) do pch_i
        pch = RampedVariable(pch_i)
        # First, must find the Tsh which gives the desired Tpr as maximum
        nlf = T -> begin
            Tsh = new_Tsh_ramp(po.Tsh, T*u"K")
            sol = solve(ODEProblem(setproperties(po, (;Tsh, pch))), odealg_chunk2)
            return ustrip(u"K", Tpr) - maximum(sol[2,:])
        end
        T_guess = ustrip(u"K", Tpr + 5.0u"K")
        # Use a secant method, provide two guess for shelf temp: 1 and 15 degrees above Tpr
        T_sh_sol = find_zero(nlf, (T_guess-4, T_guess+10),  Order1(), maxiters=10) *u"K"
        # Construct a new solution with that Tsh
        Tsh = new_Tsh_ramp(po.Tsh, T_sh_sol)
        po_i = setproperties(po, (;Tsh, pch))
        sol = solve(ODEProblem(po_i), odealg_chunk2)
        # Compute relevant quantities from the solution
        summary = summary_md_Q(sol)
        max_md = maximum(summary.md)
        dry_time = summary.t[end]
        del_t = diff(summary.t)
        weight_t = (vcat(0.0u"hr", del_t) + vcat(del_t, 0.0u"hr"))/2
        total_m = sum(summary.md .* weight_t)
        m = po.hf0*po.Ap * (po_i.ρsolution - po_i.csolid)
        ave_md = total_m/dry_time
        end_md = summary.md[end]
        max_Tf = maximum(sol[2,:])*u"K"
        return (;po=po_i, pch=pch_i, Tsh=T_sh_sol, max_md, ave_md, end_md, max_Tf, dry_time)
    end)
end

function ds_eqcap_min(func, pch_range, n_vials)
    return Table(;md=(func.(pch_range) ./ n_vials) .|> u"g/hr")
end