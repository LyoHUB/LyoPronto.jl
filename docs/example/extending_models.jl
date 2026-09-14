# # Extending LyoPronto with New Models

# LyoPronto provides a flexible framework for lyophilization modeling. This document explains 
# how to implement a new model by walking through an extension of the Pikal model with 
# microwave heating as an example.

# The core pattern for adding a new model consists of five steps:

# 1. **Define a parameter container** (`ParamObj` subtype)
# 2. **Implement the ODE right-hand-side function**
# 3. **Provide an `ODEProblem` constructor** for your parameter type
# 4. **(Optional) Add `TransformVariables` transforms** for parameter fitting
# 5. **(Optional) Add plot recipes** for visualization

# Each step leverages Julia's multiple dispatch, so the fitting machinery in 
# [`gen_sol_pd`](@ref), [`obj_pd`](@ref), [`gen_nsol_pd`](@ref), and [`objn_pd`](@ref) 
# will work with your new model automatically—provided you follow the conventions below.

# ---

# # Imports

# The following imports are used throughout this document.

using LyoPronto
using TransformVariables
using OptimizationOptimJL
using LineSearches
using ADTypes: AutoForwardDiff
using Plots
using Accessors
using LinearAlgebra: Diagonal
using ConcreteStructs: @concrete
using RecipesBase

# ---

# # 0. Write the model equations

# Before you implement a new model, you should have a clear idea of what it is.
# For this example, we will be using the following modification of the conventional Pikal 
# model for sublimation, where we add a "magic heating" term $Q_\mathrm{magic}$ that is 
# added directly to our energy balance--"magic" because it neglects all of the subtleties
# that would need to be considered for any realistic heating source.

# Heat transfer from shelf to bottom to sublimation front:
# ```math
# \begin{aligned}
# Q_\mathrm{shf} &= K_v (T_\mathrm{sh} - T_f) \\
# T_\mathrm{sub} &= T_f - \frac{Q_\mathrm{shf}}{k_\mathrm{ice}} h_\mathrm{f} 
# \end{aligned}
# ```
# Mass transfer:
# ```math
# \dot{m} = \frac{A_p}{R_p} (p_\mathrm{sub}(T_\mathrm{sub}) - p_\mathrm{ch}) 
# ```
# Overall pseudosteady energy balance and differential equation for drying progress:
# ```math
# \begin{aligned}
# \frac{d h_f}{dt} &= \frac{\dot{m}}{A_p (\rho_\mathrm{solution} - c_\mathrm{solids})} \\
# 0 &= Q_\mathrm{shf} + Q_\mathrm{magic} - \dot{m} \Delta H_\mathrm{sub}
# \end{aligned}
# ```

# In the end we have one differential equation (for ``h_f(t)``) and one more degree of 
# freedom which is fixed by our energy balance (which is an algebraic equation), so together
# we have a differential-algebraic equation (DAE) system.

# --

# # 1. Define a parameter container

# Every model in LyoPronto has a parameter container type that extends the abstract type 
# `ParamObj`. This container holds all physical, geometric, and cycle parameters needed 
# by the ODE right-hand-side function.

# ### Conventions

# - Subtype `LyoPronto.ParamObj` so the fitting machinery recognizes it.
# - Use the following names for common parameters, so that transforms used in fitting can 
#   apply across models more easily:
#   - `Rp` for the product mass transfer resistance
#   - `Kshf` for the shelf-to-frozen-product heat transfer coefficient, through a vial. In literature
#     this is commonly denoted ``K_v``, with ``v`` for vial, but in the microwave model this
#     is distinct from the vial-wall-to-frozen-product coefficient ``K_\mathrm{vw-f}``.
#   - Other names, like `pch` for ``p_\mathrm{ch}`` and `Tsh` for ``T_\mathrm{sh}``, are 
#     given without underscores for concision. Fitting for these values is less common so this matters less.
# - All the quantities which might vary from one experiment to another are captured in this 
#   struct, such as fill volume, but physical constants like the molecular weight of water
#   or enthalpy of sublimation from ice to vapor are provided separately (many by LyoPronto
#   itself, see [Physical Properties](@ref) for a listing) or use Unitful's 
#   [listing](https://juliaphysics.github.io/Unitful.jl/stable/defaultunits/#Physical-constants)
#   (e.g. `u"R"` for the gas constant).

# The struct should be declared with `@concrete terse` (from `ConcreteStructs`). The 
# `@concrete` macro ensures that all fields are concretely typed by inferring the most 
# specific type from the values provided at construction time, rather than leaving them 
# as `Any`. This is critical for performance: it avoids type instability throughout the 
# ODE right-hand-side function and the fitting pipeline, where parameters are accessed 
# repeatedly. The `terse` keyword avoids printing all the type parameters when structs are 
# shown in the REPL.

# ### Example: `ParamObjPikalMagic`

# Here is a struct with all the parameters we need:

@concrete terse struct ParamObjPikalMagic <: LyoPronto.ParamObj
    Rp
    hf0
    csolid
    ρsolution
    Kshf
    Av
    Ap
    pch
    Tsh
    Q_magic
end

# It can be helpful to define and document a helper constructor that makes sure users 
# put these parameters in the correct order, like the following which does a little 
# validation:

function ParamObjPikalMagic(tuple_of_tuples::Tuple) 
    length(tuple_of_tuples[1]) == 4 && error("Wrong tuple order given to constructor")
    length(tuple_of_tuples[2]) == 3 && error("Wrong tuple order given to constructor")
    length(tuple_of_tuples[3]) == 3 && error("Wrong tuple order given to constructor")
    popm = ParamObjPikalMagic(tuple_tuples[1]...,
    tuple_of_tuples[2]...,
    tuple_of_tuples[3]...)
    popm.pch(1.0u"hr") + 0.0u"Pa" && error("pch does not return a pressure")
    popm.Tsh(1.0u"hr") + 0.0u"K" && error("Tsh does not return an absolute temperature")
    popm.Q_magic(1.0u"hr") + 0.0u"W" && error("Q_magic does not return power")
    return popm
end

# This constructor would be used as follows: 
# ```julia
# popm = ParamObjPikalMagic((
#     (Rp, hf0, csolid, ρsolution,),
#     (Kshf, Av, Ap, ),
#     (pch, Tsh, Q_magic,),
# ))
# ```


# ---

# # 2. Implement the ODE right-hand-side function

# The ODE function has the signature `func!(du, u, params, t)` where:

# - `du` is the output derivative vector (modified in-place).
# - `u` is the current state vector (unitless).
#   - The first element of `u` and `du` should be your process completion variable, e.g. 
#     still-frozen product, that goes to 0 as drying completes.
#   - The second element of `u` and `du` should be a temperature which you want compared to 
#     experiment, e.g. the bottom center temperature where a thermocouple would be placed.
#     In the Pikal model this is fixed by an algebraic equation for pseudosteady heat and mass transfer.
# - `params` is your `ParamObj` instance.
# - `t` is the current time (unitless, in hours).

# ### Unit conventions

# - **Time** `t` is unitless but represents hours. Dimensionalize inside the function: `tn = t * u"hr"`.
# - **State** `u` is unitless but has implicit units. Assign them internally: e.g. `Tf = u[2] * u"K"`.
# - **Derivatives** `du` should be stripped back to unitless values matching the implicit units: 
#   `du[2] = ustrip(u"K/hr", dTf)`.

# ### Example: `pikal_mw!`

# The extended Pikal model has two state variables: `[hf, Tf]` (remaining frozen layer 
# thickness, product temperature). It is implemented as a DAE (differential-algebraic 
# equation) with mass matrix `Diagonal([1.0, 0.0])`, where the algebraic constraint 
# enforces the pseudosteady-state energy balance.

# Following the Pikal model pattern, the physics is split into two functions:

# 1. A method for the **helper function** [`calc_md_Q`](@ref), dispatched on the new 
#    `ParamObjPikalMagic` type, that computes physical quantities 
#    (mass flow, heat transfer terms) and returns them as a named tuple. This function 
#    can be reused independently for diagnostics or post-solution analysis.
# 2. The **RHS function** (`pikal_mw!`) that calls the helper, then writes derivatives 
#    and algebraic residuals to `du`.

# #### Helper: `calc_md_Q`

# This function unpacks parameters, dimensionalizes state and time, computes all heat 
# and mass transfer terms, and returns a named tuple.

# It is *very important* that this be dispatched on the new parameter struct (`ParamObjPikalMagic`
# in this case) to ensure it is not mixed up with the function as defined for other sets of 
# model equations.

# This function should return at least the mass flow rate as `md` and the shelf-to-frozen-product
# heat transfer as `Q_shf`, in order to interoperate with fitting and plotting functions.
# Feel free to include as many other computed quantities as will be of use (and would be 
# annoying to compute): for example, ``T_\mathrm{sub}`` can be computed from known `Q_shf` 
# and `Tf`, but if you are interested in plotting it, you can keep your plotting logic in 
# sync with the model logic by returning `Tsub` in the named tuple here.

@inline function calc_md_Q(u, po::ParamObjPikalMagic, t)
    (;Rp, hf0, Kshf, Av, Ap, pch, Tsh) = po
    
    td = t * u"hr"
    hf = u[1] * u"cm"
    Tf = u[2] * u"K"
    hd = hf0 - hf
    ## Shelf heat transfer
    Q_shf = Kshf(pch(td)) * Av * (Tsh(td) - Tf) |> u"W"
    ## Microwave heating (NEW term)
    Q_magic = po.Q_magic(t) |> u"W"
    ## Sublimation mass transfer
    Tsub = Tf - Q_shf / k_ice / Ap * hf
    delta_p = calc_psub(Tsub) - pch(td)
    md = -Ap * delta_p / Rp(hd) |> u"g/hr"
    return (; md, Q_shf, Q_mw, Tsub) # Only `md` and `Q_shf` are crucial to the RHS function
end

# #### RHS: `pikal_mw!`

# The RHS function calls the helper, then computes and writes the derivative and residuals.

function pikal_mw!(du, u, params, t)
    (;md, Q_shf, Q_mw) = calc_md_Q_mw(u, params, t)
    (; csolid, ρsolution, Ap) = params

    Q_sub = uconvert(u"W", md * ΔHsub)
    dhf_dt = min(0.0u"cm/hr", md / (ρsolution - csolid) / Ap |> u"cm/hr")

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Q_sub + Q_shf + Q_mw)
end

# Because this modeal is a DAE system, wrap the RHS in an `ODEFunction` with a `Diagonal` mass matrix, with entries
# of `1.0` for differential equations and `0.0` for algebraic equations. 

const pikal_mw_f = ODEFunction(pikal_mw!, mass_matrix=Diagonal([1.0, 0.0]))

# ### For your model

# If your model is a DAE (like the Pikal model), wrap the RHS in an `ODEFunction` with 
# a `mass_matrix` argument, then pass the `ODEFunction` to `ODEProblem`. If it is a set of
# ODEs without any algebraic constraints, the RHS function as defined can be passed directly
# to the `ODEProblem` constructor.

# ---

# # 3. Define an `ODEProblem` constructor

# The fitting functions [`gen_sol_pd`](@ref) and [`gen_nsol_pd`](@ref) call 
# `ODEProblem(param_obj)` to construct the problem. By defining a method for your 
# `ParamObj` subtype, your model can be integrated with the entire fitting pipeline.

# ### Required helper methods

# Before defining `ODEProblem`, you need two helpers:

# #### `calc_u0(po::YourParamObj)`

# Returns the initial condition vector `u0` as a unitless `Vector{Float64}`. 
# This needs to be kept in sync with the RHS function defined previously.
# Note also that the default DAE initialization used to satisfy algebraic constraints,
# `BrownFullBasicInit()`, will treat the initial conditions for algebraic variables as 
# a guess, so algebraic variables need not be exact for this function.

function calc_u0(po::ParamObjPikalMagic)
    return [ustrip(u"cm", po.hf0), ustrip(u"K", float(po.Tsh(0u"s")))]
end

# #### `get_tstops(po::YourParamObj)`

# This is used to compute a sorted, unique vector of time points where cycle conditions 
# aren't differentiable (e.g. at the corner where ramp goes to flat setpoint hold). 
# These are used as `tstops` in the ODE solver for accuracy.

# The method of this function which accepts a `Tuple` will automatically call [`LyoPronto.extract_ts`](@ref)
# which identifies all the corners of a [`RampedVariable`](@ref) and all the time points of a 
# [`LinearInterpolation`](https://docs.sciml.ai/DataInterpolations/stable/methods/#Linear-Interpolation).

function get_tstops(po::ParamObjPikalMagic)
    get_tstops((po.Tsh, po.pch, po.Q_magic))
end


# #### `ODEProblem(po::YourParamObj; u0=calc_u0(po), tspan=(0.0, 1000.0))`

# In the type parameters of 
# [`ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction),
# the `true` informs the solver that the derivatives are computed in-place in a passed vector,
# and `FullSpecialize` ensures that compilation is fully specialized on all types. 
# Using `FullSpecialize` incurs a greater compilation cost to speed up calculation, but more 
# importantly it avoids some kinds of problems with automatic differentation (particularly
# [this issue](https://github.com/SciML/OrdinaryDiffEq.jl/issues/3381), at time of writing).

function ODEProblem(po::ParamObjPikalMagic; u0=calc_u0(po), tspan=(0.0, 1000.0))
    tstops = get_tstops(po)
    return ODEProblem{true, SciMLBase.FullSpecialize}(
        pikal_mw_f, u0, tspan, po;
        tstops=tstops, callback=end_drying_callback,
        initializealg=BrownFullBasicInit(), dt=0.1
    )
end

# Key points:
# - The first argument to `ODEProblem` is your RHS `ODEFunction` (or bare function for explicit ODEs).
# - `callback=end_drying_callback` terminates integration when drying is complete 
#   (frozen layer thickness approaches zero, i.e. ≈ 1e-10).
# - For DAEs, `initializealg=BrownFullBasicInit()` provides consistent initial conditions.
# - `tspan` should be very large, e.g. 1000 hours—the callback will stop simulation when drying is done.
# - Use [`get_tstops`](@ref) as defined above is used to ensure simulation treats all the  
#   non-smooth points.
# - `dt=0.1` is a conservative estimate; time steps often exceed hours.
# - This is a natural place to add any other 
#   [keyword arguments](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/) 
#   for the solver, such as tolerances.

# ---

# # 4. (Optional) Add `TransformVariables` transforms for fitting

# LyoPronto uses the `TransformVariables` package to map unconstrained optimization 
# parameters to physically meaningful ranges. Convenience functions already exist for 
# common parameters:

# - [`K_transform_basic`](@ref) — transforms for `Kshf` (shelf heat transfer coefficient)
# - [`Rp_transform_basic`](@ref) — transforms for `Rp` (product resistance)

# ### How transforms work

# A transform maps a flat `Vector{Float64}` to a `NamedTuple` of parameter values:

trans_K = K_transform_basic(5.0u"W/m^2/K")
# Maps a scalar to (; Kshf = ConstPhysProp(5.0u"W/m^2/K"))

# To add a transform for a new parameter, compose 
# [scalar transforms](https://www.tamaspapp.eu/TransformVariables.jl/stable/#Scalar-transforms)
# with appropriate units for each parameter, and nest them inside a transform to named tuple.

# ### How fitting uses transforms

# The fitting functions expect a tuple `tpf = (transform, param_objs, fitdats)`. 
# Internally they:

# 1. Call `transform(tr, fitlog)` to get a `NamedTuple` of parameters.
# 2. Call `setproperties(po, fitprm)` to merge fitted params into the base `ParamObj`.
# 3. Call `ODEProblem(new_po)` to construct the ODE, then solve the ODE.
# 4. Compare the solution to data in a `PrimaryDryFit` .

# Because `setproperties` (from `ConstructionBase`) works on any struct, and `ODEProblem` 
# dispatches on your `ParamObj` subtype, **no additional code is needed** for fitting to work.

# ---

# # 5. (Optional) Add plot recipes

# LyoPronto uses `RecipesBase` to provide `plot` recipes for custom types. 

# The plot recipe [`modconvtplot`](@ref) will plot the 2nd variable of your state vector `u`
# (defined for the ODE RHS function) as a temperature. If you would like to plot other 
# quantities specific to your model, you can add [recipes](https://docs.juliaplots.org/latest/RecipesBase/types/#User-Recipes-2) like the following:

@userplot PikalMagicQPlot
@recipe function f(pmqp::PikalMagicQPlot)
    sol = pqmp.args
    summary = summary_md_Q(sol)
    @series begin
        color --> "orange"
        return summary.t, summary.Q_shf
    end
    @series begin
        color --> "dark orange"
        return summary.t, summary.Q_magic
    end
end

# For a solution `sol` with the above model, this would then be called as
# ```julia
# pikalmagicqplot(sol)
# ```

# ---

# # Summary checklist

# To add a new model to LyoPronto, implement:

# | Step | Required? | Purpose |
# |------|-----------|---------|
# | `ParamObj` subtype | Yes | Container for all model parameters |
# | ODE RHS function | Yes | Model physics |
# | `ODEFunction` with mass matrix | Yes if algebraic equations | For DAE models (e.g., Pikal-style) |
# | `calc_u0` method | Yes | Initial conditions |
# | `get_tstops` method | Yes | Solver accuracy at ramp transitions |
# | `ODEProblem` method | Yes | Problem construction for solving and fitting |
# | `TransformVariables` transforms | Optional | Parameter fitting support |
# | Plot recipes | Optional | Visualization |

# Once these are in place, the existing functions 
# [`gen_sol_pd`](@ref), [`obj_pd`](@ref), [`gen_nsol_pd`](@ref), [`objn_pd`](@ref), 
# [`nls_pd`](@ref), and [`nls_pd!`](@ref) will work with your model automatically.
